/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/


#include "cuda_runtime.h"
#include <nvToolsExt.h>

#include "as/GPUWorker.h"
#include "asUtils/macros.h"
#include "asUtils/timer.h"
#include "as/ServerManagement.h"
#include "as/asBuffer.h"
#include "as/asTypes.h"
#include "asCore/asCoreHelpers.h"
#include "asCore/Transformer.h"
#include "asCore/Interpolator.h"
#include "as/RingArray.h"
#include "as/WorkItem.h"
#include "asUtils/Logger.h"

#include "config.h"

#include <iostream>
#include <cassert>

using namespace Log;

/* Constructor */
as::GPUWorker::GPUWorker(ServerManagement& S_mngt,
		int deviceId, unsigned atomBufferSize, unsigned dofBufferSize) :
		Worker(S_mngt), _deviceId(deviceId),
		_atomBufferSize(atomBufferSize), _dofBufferSize(dofBufferSize) {}

/* Destructor */
as::GPUWorker::~GPUWorker ()
{
	global_log->info() << std::setw(5) << " " << "Delete GPUWorker of Device: " << _deviceId;
	if(_queue.sleeps() == true || _queue.terminates() == false) {
		*global_log << " Error" << std::endl;
	} else {
		*global_log << " Ok" << std::endl;
	}
};

/****************************
 * public member functions
 ****************************/

void as::GPUWorker::run () {

	/* Set the device to work with */
	CUDA_CHECK(cudaSetDevice(_deviceId));

#ifndef KERNELTIMING
	/* Initialize Pipeline */
	static const unsigned numStages = 5;
	bool predicates[2][numStages];
	for (unsigned i = 0; i<numStages; ++i) {
		predicates[0][i] = false;
		predicates[1][i] = false;
	}
	RingArray<WorkerItem*> stagesMngt(numStages);
	RingArray<Protein*> LigMngt(numStages);
	unsigned pipeIdx[2] = {0,1};
	int numItemsInPipe = 0;


	/* Streams */
	cudaStream_t streams[4];
	for (int i = 0; i<4; ++i) {
		CUDA_CHECK(cudaStreamCreate(&streams[i]));
	}

	/* Events */
	cudaEvent_t events[7];
	for (int i = 0; i<7; ++i) {
		CUDA_CHECK(cudaEventCreate(&events[i], cudaEventDisableTiming));
	}

	/* Initialize Input/Output Buffers */
	Comp1_HD<DOF, DEVONLY> d_dof0(_dofBufferSize);
	d_dof0.initDevice();
	Comp1_HD<DOF, DEVONLY> d_dof1(_dofBufferSize);
	d_dof1.initDevice();
	Comp1_HD<DOF, DEVONLY>* d_dof[2];
	d_dof[0] = &d_dof0;
	d_dof[1] = &d_dof1;

	Comp3_HD<float, DEVONLY> d_trafoLig(_atomBufferSize);
	d_trafoLig.initDevice();

	Comp5_HD<float, DEVONLY> d_potLig0(_atomBufferSize);
	d_potLig0.initDevice();
	Comp5_HD<float, DEVONLY> d_potLig1(_atomBufferSize);
	d_potLig1.initDevice();
	Comp5_HD<float, DEVONLY>* d_potLig[2];
	d_potLig[0] = &d_potLig0;
	d_potLig[1] = &d_potLig1;

	//Debug
//	Comp3_HD<float, HOST> h_trafoLig(_atomBufferSize);
//	h_trafoLig.initHost();
	//DEBUG
//	Comp5_HD<float, HOST> h_potLig(_atomBufferSize);
//	h_potLig.initHost();

	Comp1_HD<float, HOST_PINNED> hd_res0(14*_dofBufferSize);
	hd_res0.initHost();
	hd_res0.initDevice();
	Comp1_HD<float, HOST_PINNED> hd_res1(14*_dofBufferSize);
	hd_res1.initHost();
	hd_res1.initDevice();
	Comp1_HD<float, HOST_PINNED>* hd_res[2];
	hd_res[0] = &hd_res0;
	hd_res[1] = &hd_res1;

	/* Create local Transformer and Interpolator object */
	asCore::Transformer TF;
	asCore::Interpolator IP;

	/* Create SimParam obj */
	/* ToDo: object of this type should be added to the dataMngt in global and per-client manner */
//	SimParam simPar;
//	simPar.dielec = variable;				/** type of dielectric constant */
//	simPar.epsilon = 15;					/** dielectric constant */
//	simPar.ffelec = FELEC/simPar.epsilon;	/** precomputed factor felec/epsilon */
////	simPar.useSwi = false;					/** using switching potential */
////	simPar.swiOn = 0;						/** min. switching potential distance */
////	simPar.swiOff= 0;						/** max. switching potential distance */
//	simPar.useRecGrad = false;				/** using Receptor gradients */
//	simPar.usePot = true;					/** use Potential grid */

	//DEBUG
//	int count = 0;


	while (true) {

		/* reset the predicates for the actual iteration*/
		for (unsigned i = 0; i<numStages; ++i) {
			predicates[pipeIdx[0]][i] = false;
		}

		/* get a workItem from the queue */
		if (numItemsInPipe == 0) {
			WorkerItem* item;
//			nvtxRangePushA("GPUWorker BlockedRemove");
			item = _queue.removeCondBlocked();
//			nvtxRangePop();
			if (item == NULL) {
				/* the worker was signaled to terminate */
				break;
			}
			assert(item->size() > 0);
			/* get device local ids */
			item->setDevLocGridId(deviceLocalGridID(item->globGridId()));
			item->setDevLocLigId(deviceLocalProteinID(item->globLigId()));
			item->setDevLocRecId(deviceLocalProteinID(item->globRecId()));

			++numItemsInPipe;
			/* signal that stage 0 is to executed within the current iteration */
			predicates[pipeIdx[0]][0] = true;
			stagesMngt.push(item);
			LigMngt.push(getProtein(item->globLigId()));
		} else {
			WorkerItem* item;
//			nvtxRangePushA("GPUWorker Remove");
			item = _queue.removeUnblocked();
//			nvtxRangePop();
			if (item != NULL) {
				assert(item->size() > 0);
				/* get device local ids */
				item->setDevLocGridId(deviceLocalGridID(item->globGridId()));
				item->setDevLocLigId(deviceLocalProteinID(item->globLigId()));
				item->setDevLocRecId(deviceLocalProteinID(item->globRecId()));

				++numItemsInPipe;
				/* signal that stage 0 is to executed within the current iteration */
				predicates[pipeIdx[0]][0] = true;
				stagesMngt.push(item);
				LigMngt.push(getProtein(item->globLigId()));
			} else {
				stagesMngt.rotate();
				LigMngt.rotate();
			}
		}

//		std::cout << "begin" << std::endl;
//		for (int j = 0; j < 2; ++j) {
//			for (unsigned i = 0; i<numStages; ++i) {
//				std::cout << predicates[pipeIdx[j]][i] ;
//			}
//			std::cout << std::endl;
//		}

		/* stage 0
		 * Copy data to device */

		if (predicates[pipeIdx[0]][0] == true)
		{
			const static unsigned stageId = 0;
			const WorkerItem* it = stagesMngt.get(stageId);
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;

			cudaVerify(cudaStreamWaitEvent(streams[0], events[2], 0));
			cudaVerify(cudaMemcpyAsync(d_dof[pipeIdx[0]]->d_data(), it->DOFBuffer(),
					it->size()*sizeof(DOF),	cudaMemcpyHostToDevice, streams[0]));
			cudaVerify(cudaEventRecord(events[0], streams[0]));

		}

		/* stage 1
		 * Transform coordinates and perform interpolation and neigborlist calculations */

		/* check if stage 0 was executed in last iteration */
		if (predicates[pipeIdx[1]][0] == true)
		{
			const static unsigned stageId = 1;
			const WorkerItem* it = stagesMngt.get(stageId);
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;

			/* Device: Wait for completion of copyH2D of DOFs to complete */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[0], 0));

			const unsigned numEl = it->size()*LigMngt.get(stageId)->nAtoms();
//			std::cout << "GPUWorker.cpp: " << "numEl " << numEl << " _atomBufferSize " << _atomBufferSize << std::endl;
			assert(numEl <= _atomBufferSize);

			/* Perform cuda kernel calls */
			TF.calcAndSetGridSize(numEl);
			TF.d_DOF2Pos(it->devLocLigId(),it->size(),d_dof[pipeIdx[1]], &d_trafoLig, streams[2]);

			//DEBUG
//			if(count++ == 18494) {
//				h_trafoLig.set_d_x(d_trafoLig.d_x());
//				h_trafoLig.set_d_y(d_trafoLig.d_y());
//				h_trafoLig.set_d_z(d_trafoLig.d_z());
//				cudaDeviceSynchronize();
//				h_trafoLig.cpyD2H();
//				cudaDeviceSynchronize();
//				h_trafoLig.printEl(numEl, 0);
////					return 0;
//			}

			/* Device: Signal event when transformation has completed */
			cudaVerify(cudaEventRecord(events[2], streams[2]));
			/* Device: Wait for completion of reduction of the previous round */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[5+pipeIdx[1]], 0));

			/* Perform cuda kernel calls */
			IP.calcAndSetGridSize(numEl);
			IP.d_PotForce<asCore::built_in>(it->devLocGridId(), it->devLocLigId(), it->size(), &d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);
//			IP.d_PotForce<manual>(it->devLocGridId(), it->devLocLigId(), it->size(), &d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);

			//DEBUG
//			if(count++ == 18494) {
//				IP.d_PotForce<manual>(it->devLocGridId(), it->devLocLigId(), it->size(), &d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);
//				h_potLig.set_d_x(d_potLig[pipeIdx[1]]->d_x());
//				h_potLig.set_d_y(d_potLig[pipeIdx[1]]->d_y());
//				h_potLig.set_d_z(d_potLig[pipeIdx[1]]->d_z());
//				h_potLig.set_d_v(d_potLig[pipeIdx[1]]->d_v());
//				h_potLig.set_d_w(d_potLig[pipeIdx[1]]->d_w());
//				h_potLig.cpyD2H();
//				cudaDeviceSynchronize();
////				h_potLig.printEl(numEl, 0);
////					return 0;
//			}

			IP.d_NLPotForce<false>(it->devLocGridId(), it->devLocRecId(), it->devLocLigId(),it->size(),
					&d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);

			//DEBUG
//			if(count++ == 18494) {
//				IP.d_NLPotForce<false>(it->devLocGridId(), it->devLocRecId(), it->devLocLigId(),it->size(),
//					&d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);
//				h_potLig.set_d_x(d_potLig[pipeIdx[1]]->d_x());
//				h_potLig.set_d_y(d_potLig[pipeIdx[1]]->d_y());
//				h_potLig.set_d_z(d_potLig[pipeIdx[1]]->d_z());
//				h_potLig.set_d_v(d_potLig[pipeIdx[1]]->d_v());
//				h_potLig.set_d_w(d_potLig[pipeIdx[1]]->d_w());
//				h_potLig.cpyD2H();
//				cudaDeviceSynchronize();
//				h_potLig.printEl(numEl, 0);
////					return 0;
//			}

			/* Device: Signal event when force and energy calc. has completed */
			cudaVerify(cudaEventRecord(events[3+pipeIdx[1]], streams[2]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][1] = true;
		}

		/* stage 2
		 * Reduce the result of PotForce calc. to force gradients */

		/* check if stage 1 was executed in last iteration */
		if (predicates[pipeIdx[1]][1] == true)
		{
			const static unsigned stageId = 2;
			const WorkerItem* it = stagesMngt.get(stageId);

			const unsigned nAtoms = LigMngt.get(stageId)->nAtoms();
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;

			/* Device: Wait for completion of PotForce calc. to complete */
			cudaVerify(cudaStreamWaitEvent(streams[3], events[3+pipeIdx[0]], 0));

			/* Perform cuda kernel calls */
			/* old version */
//			TF.d_partForce2Grad(it->devLocLigId(), it->size(), d_potLig[pipeIdx[0]], hd_res[pipeIdx[0]], streams[3]);

			/* new version */
			TF.d_partForce2GradAll(it->devLocLigId(), it->size(), nAtoms, d_potLig[pipeIdx[0]], hd_res[pipeIdx[0]], streams[3]);

			/* Device: Signal event when reduction has completed */
			cudaVerify(cudaEventRecord(events[5+pipeIdx[0]], streams[3]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][2] = true;
		}

		/* stage 3
		 * copy partial results back from the device */

		/* check if stage 2 was executed in last iteration */
		if (predicates[pipeIdx[1]][2] == true)
		{
			const static unsigned stageId = 3;
			const WorkerItem* it = stagesMngt.get(stageId);
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;

			/* Device: Wait for completion of reduction calc. to complete */
			cudaVerify(cudaStreamWaitEvent(streams[1], events[5 + pipeIdx[1]], 0));
			cudaVerify(cudaMemcpyAsync(hd_res[pipeIdx[1]]->h_data(), hd_res[pipeIdx[1]]->d_data(),
					14*it->size()*sizeof(float), cudaMemcpyDeviceToHost, streams[1]));

			/* Device: Signal event when transfere has completed */
			cudaVerify(cudaEventRecord(events[1], streams[1]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][3] = true;
		}

//		/* stage 5
//		 * This stage is chronologically called after stage 4 because it waits for completion of
//		 * stage 4 belonging the previous iteration.
//		 * However, it has to be called before stage 4 in order process data before the host thread
//		 * gets synchronized to wait for completion of device tasks that were scheduled in previous stages.
//		 *
//		 * This stage is responsible for copying results to the EnGrad buffer of the respective work item.*/
//
//		/* check if stage 4 was executed in last iteration */
//		if (predicates[pipeIdx[1]][4] == true)
//		{
//			const static unsigned stageId = 5;
//			const WorkerItem* it = stagesMngt.get(stageId);
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;
//
//
//
//			/* signal that this stage was executed within the current iteration */
//			/* not needed because not stage waits for completion but good for monitoring
//			 * what is happening*/
//
//
//
//			predicates[pipeIdx[0]][5] = true;
//			/* signal that one item has been passed the last stage */
//			--numItemsInPipe;
//		}

		/* stage 4
		 * Performs final reduction and copies result directly to EnGrad buffer of work item */

		/* check if stage 3 was executed in last iteration */
		if (predicates[pipeIdx[1]][3] == true)
		{
			const static unsigned stageId = 4;
			WorkerItem* const it = stagesMngt.get(stageId);
//			std::cout << "Stage " << stageId << " processing" << std::endl;
//			std::cout << "NumEl2Process " << it->size() << std::endl;

			/* Host: Wait for completion of data transfer to complete */
			cudaVerify(cudaEventSynchronize(events[1]));
			nvtxRangePushA("Host");
			TF.h_finalForce2Grad(LigMngt.get(stageId),it->size(), it->DOFBuffer(), hd_res[pipeIdx[0]], it->EnGradBuffer());
			nvtxRangePop();
			/* Signal that result is in buffer */
			it->setReady();

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][4] = true;

			/* signal that one item has been passed the last stage */
			--numItemsInPipe;
		}

		assert(numItemsInPipe >= 0);

//		std::cout << numItemsInPipe << "items in pipe" << std::endl;
		std::swap(pipeIdx[0], pipeIdx[1]);

//		std::cout << std::endl;
//		std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	}

	/* Free/Destroy resources */
	for (int i = 0; i<4; ++i) {
		CUDA_CHECK(cudaStreamDestroy(streams[i]));
	}
	for (int i = 0; i<7; ++i) {
		CUDA_CHECK(cudaEventDestroy(events[i]));
	}

#endif

#ifdef KERNELTIMING

	using namespace std;
	const int numSteps = 80;
	const int powMax = 4;
	const int powMin = 0;
	float incr = float(powMax - powMin)/(numSteps-1);
	vector<float> exponents(numSteps);
	set<unsigned> chunks;

	for(int i = 0; i < numSteps-1; ++i) {
		exponents[i] = float(powMin) + i*incr;
//		cout << exponents[i] << endl;
	}
	exponents[numSteps-1] = powMax;


	for(int i = 0; i < numSteps; ++i) {
		chunks.insert(unsigned(pow(10.0, exponents[i])));
	}

	int precisionSetting = cerr.precision( );
	ios::fmtflags flagSettings = cerr.flags();
	cerr.setf(ios::scientific | ios::showpos);
	cerr.precision(5);

	unsigned pipeIdx[2] = {0,1};
	cudaEvent_t start, stop;
	cudaVerify(cudaEventCreate(&start));
	cudaVerify(cudaEventCreate(&stop));

	/* Create local Transformer and Interpolator object */
	Transformer TF;
	Interpolator IP;

	/* Create SimParam obj */
	/* ToDo: object of this type should be added to the dataMngt in global and per-client manner */
//	SimParam simPar;
//	simPar.dielec = variable;				/** type of dielectric constant */
//	simPar.epsilon = 15;					/** dielectric constant */
//	simPar.ffelec = FELEC/simPar.epsilon;	/** precomputed factor felec/epsilon */
////	simPar.useSwi = false;					/** using switching potential */
////	simPar.swiOn = 0;						/** min. switching potential distance */
////	simPar.swiOff= 0;						/** max. switching potential distance */
//	simPar.useRecGrad = false;				/** using Receptor gradients */
//	simPar.usePot = true;					/** use Potential grid */



	/* Only one item available */
	/* get a workItem from the queue */

	WorkerItem* item;
	nvtxRangePushA("GPUWorker BlockedRemove");
	item = _queue.removeCondBlocked();
	nvtxRangePop();

	assert(item->size() > 0);
	/* get device local ids */
	item->setDevLocGridId(deviceLocalGridID(item->globGridId()));
	item->setDevLocLigId(deviceLocalProteinID(item->globLigId()));
	item->setDevLocRecId(deviceLocalProteinID(item->globRecId()));

	const Protein* lig = getProtein(item->globLigId());


	for (auto chunk : chunks) {
		double tDOF2Coord, tPotForce_m, tPotForce_bi, tNL, t_h_Force2Grad, t_d_Force2Grad, tCpyH2D, tCpyD2H, t5CompD2H, t3CompH2D;
		tDOF2Coord = tPotForce_m = tPotForce_bi = tNL = t_h_Force2Grad  = t_d_Force2Grad = tCpyH2D = tCpyD2H = t5CompD2H = t3CompH2D = 0.0;

		unsigned atomBufferSize = chunk * lig->ntotAtoms();
		/* Initialize Input/Output Buffers */
		Comp1_HD<DOF, DEVONLY> d_dof(chunk);
		d_dof.initDevice();

		Comp3_HD<float, DEVONLY> d_trafoLig(atomBufferSize);
		d_trafoLig.initDevice();

		Comp5_HD<float, DEVONLY> d_potLig(atomBufferSize);
		d_potLig.initDevice();

		Comp1_HD<float, HOST_PINNED> hd_res(14*atomBufferSize);
		hd_res.initHost();
		hd_res.initDevice();

		Comp3_HD<float, HOST_PINNED> h_trafoLig(atomBufferSize);
		h_trafoLig.initHost();
		h_trafoLig.set_d_x(d_trafoLig.d_x());
		h_trafoLig.set_d_y(d_trafoLig.d_y());
		h_trafoLig.set_d_z(d_trafoLig.d_z());

		Comp5_HD<float, HOST_PINNED> h_potLig(atomBufferSize);
		h_potLig.initHost();
		h_potLig.set_d_x(d_potLig.d_x());
		h_potLig.set_d_y(d_potLig.d_y());
		h_potLig.set_d_z(d_potLig.d_z());
		h_potLig.set_d_v(d_potLig.d_v());
		h_potLig.set_d_w(d_potLig.d_w());



		unsigned numStruc = item->size();
		uint numRounds = (numStruc + chunk - 1) / chunk;
		for (uint rnd = 0, dofIdx = 0; rnd < numRounds; ++rnd, dofIdx+=chunk) {
			int chunkSize = (numStruc - dofIdx) > chunk ? chunk : (numStruc - dofIdx);
			/* stage 0
			 * Copy data to device */
			if (true)
			{
//				if (chunkSize == 5)
//					cout << "DOFidx: " << dofIdx << endl;
				cudaInitTimer(start);
				cudaVerify(cudaMemcpyAsync(d_dof.d_data(), item->DOFBuffer() + dofIdx,
						chunkSize*sizeof(DOF),	cudaMemcpyHostToDevice));
				tCpyH2D += cudaGetTimer(start, stop);
			}

			/* stage 1
			 * Transform coordinates and perform interpolation and neigborlist calculations */

			/* check if stage 0 was executed in last iteration */
			if (true)
			{

				const unsigned numEl = chunkSize*lig->nAtoms();
				assert(numEl <= atomBufferSize);

				/* Perform cuda kernel calls */
				TF.calcAndSetGridSize(numEl);

				cudaInitTimer(start);
				TF.d_DOF2Pos(item->devLocLigId(),chunkSize, &d_dof, &d_trafoLig);
				tDOF2Coord += cudaGetTimer(start, stop);

				h_trafoLig.cpyD2H();
				cudaInitTimer(start);
				h_trafoLig.cpyH2D();
				t3CompH2D += cudaGetTimer(start, stop);

				/* Perform cuda kernel calls */
				IP.calcAndSetGridSize(numEl);

				cudaInitTimer(start);
				IP.d_PotForce<built_in>(item->devLocGridId(), item->devLocLigId(), chunkSize, &d_trafoLig, &d_potLig);
				tPotForce_bi += cudaGetTimer(start, stop);

				cudaInitTimer(start);
				IP.d_PotForce<manual>(item->devLocGridId(), item->devLocLigId(), chunkSize, &d_trafoLig, &d_potLig);
				tPotForce_m += cudaGetTimer(start, stop);

				cudaInitTimer(start);
				IP.d_NLPotForce<false>(item->devLocGridId(), item->devLocRecId(), item->devLocLigId(),chunkSize,
					&d_trafoLig, &d_potLig);
				tNL += cudaGetTimer(start, stop);

				cudaInitTimer(start);
				h_potLig.cpyH2D();
				t5CompD2H += cudaGetTimer(start, stop);

			}

			/* stage 2
			 * Reduce the result of PotForce calc. to force gradients */

			/* check if stage 1 was executed in last iteration */
			if (true)
			{
				/* Perform cuda kernel calls */
				cudaInitTimer(start);
				/* old version */
//				TF.d_partForce2Grad(item->devLocLigId(), chunkSize, &d_potLig, &hd_res);
				/* new version */
				TF.d_partForce2GradAll(item->devLocLigId(), chunkSize, lig->nAtoms(), &d_potLig, &hd_res);
				t_d_Force2Grad += cudaGetTimer(start, stop);
			}

			/* stage 3
			 * copy partial results back from the device */

			/* check if stage 2 was executed in last iteration */
			if (true)
			{
				/* Device: Wait for completion of reduction calc. to complete */
				cudaInitTimer(start);
				cudaVerify(cudaMemcpyAsync(hd_res.h_data(), hd_res.d_data(),
						14*chunkSize*sizeof(float), cudaMemcpyDeviceToHost));
				tCpyD2H += cudaGetTimer(start, stop);

			}

			/* stage 4
			 * Performs final reduction and copies result directly to EnGrad buffer of work item */

			/* check if stage 3 was executed in last iteration */
			if (true)
			{
				Timer timer;
				initTimer(&timer);
				TF.h_finalForce2Grad(lig,chunkSize, item->DOFBuffer() + dofIdx, &hd_res, item->EnGradBuffer() + dofIdx);
				t_h_Force2Grad += getTimer(&timer);
			}

			std::swap(pipeIdx[0], pipeIdx[1]);
		}

		std::cerr	<< "TIME" << "\t"
					<< chunk << "\t"
					<< tCpyH2D << "\t"
					<< tDOF2Coord << "\t"
					<< tPotForce_bi << "\t"
					<< tPotForce_m << "\t"
					<< tNL << "\t"
					<< t_d_Force2Grad << "\t"
					<< tCpyD2H << "\t"
					<< t_h_Force2Grad << "\t"
					<< t3CompH2D << "\t"
					<< t5CompD2H << std::endl;
	}

	cerr.precision(precisionSetting);
	cerr.flags(flagSettings);


	/* Signal that result is in buffer */
	item->setReady();

	/* makes the worker sleeping */
	item = _queue.removeCondBlocked();
#endif


	global_log->info() << std::setw(5) << " " << "GPUWorker " << _deviceId <<" terminated" << std::endl;

}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/

inline int as::GPUWorker::deviceLocalProteinID(const int& globalId) {
	return _S_mngt.DataMngt()->deviceLocalProteinID(_deviceId, globalId);
}

inline int as::GPUWorker::deviceLocalGridID(const int& globalId) {
	return _S_mngt.DataMngt()->deviceLocalGridID(_deviceId, globalId);
}


