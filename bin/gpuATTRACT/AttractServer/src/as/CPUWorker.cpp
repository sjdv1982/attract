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

#include "as/CPUWorker.h"
#include "asUtils/timer.h"
#include "as/ServerManagement.h"
#include "as/asBuffer.h"
#include "as/asTypes.h"
#include "asCore/Transformer.h"
#include "asCore/Interpolator.h"
#include "asUtils/Logger.h"

#include "config.h"

#include <iostream>
#include <cassert>

using Log::global_log;

/* Constructor */
as::CPUWorker::CPUWorker(ServerManagement& S_mngt,
		unsigned atomBufferSize) :
		Worker(S_mngt),
		_id(-1),
		_atomBufferSize(atomBufferSize) {}

/* Destructor */
as::CPUWorker::~CPUWorker ()
{
	global_log->info() << std::setw(5) << " " << "Delete CPUWorker: " << _id;
	if(_queue.sleeps() == true || _queue.terminates() == false) {
		*global_log << " Error" << std::endl;
	} else {
		*global_log << " Ok" << std::endl;
	}
};

/****************************
 * public member functions
 ****************************/
void as::CPUWorker::run ()
{

#ifndef KERNELTIMING
	/* Initialize Input/Output Buffers */

	Comp3_HD<float, HOSTONLY> h_trafoLig(_atomBufferSize);
	h_trafoLig.initHost();

	Comp5_HD<float, HOSTONLY> h_potLig(_atomBufferSize);
	h_potLig.initHost();

	//DEBUG
//	Comp5_HD<float, HOSTONLY> h_potLigNL(_atomBufferSize);
//	h_potLigNL.initHost();

	/* Create local Transformer and Interpolator object */
	asCore::Transformer TF;
	asCore::Interpolator IP;

	while (true) {
		/* get a workItem from the queue */
		WorkerItem* item;
		nvtxRangePushA("CPUWorker BlockedRemove");
		item = _queue.removeCondBlocked();
		nvtxRangePop();
		if (item == NULL) {
			/* the worker was signaled to terminate */
			break;
		}

		assert(item->size() > 0);

		/* each DOF is processed sequentially */

		//DEBUG
//		static int count = -1;
//		count++;

		const Protein* rec = getProtein(item->globRecId());
		const Protein* lig = getProtein(item->globLigId());
		const GridUnion* grid = getGridUnion(item->globGridId());
		const AttrParamTable* table = getParamTable();
		const SimParam* simPar = getSimParam();

		for (unsigned i = 0; i < item->size(); ++i) {
			nvtxRangePushA("CPUWorker Processing");
			const DOF& dof = item->DOFBuffer()[i];
			EnGrad& engrad = item->EnGradBuffer()[i];
			asUtils::RotMatf rotMat;

			asCore::euler2rotmat(dof.ang.x, dof.ang.y, dof.ang.z, rotMat);
			TF.h_DOF2Pos(lig, dof, rotMat, &h_trafoLig);

			//DEBUG
//			if(count++ == 18494 &&  true) {
//				h_trafoLig.printEl(lig->numAtoms(), 0);
//				return;
//			}

			IP.h_PotForce(grid,lig, &h_trafoLig, &h_potLig);

			//DEBUG
//			if(count++ == 18494 &&  true) {
//				h_potLigNL.resetHostData();
//				IP.h_PotForce(grid,lig, &h_trafoLig, &h_potLigNL);
//				h_potLigNL.printEl(lig->numAtoms(), 0);
////				return;
//			}


			IP.h_NLPotForce(grid, rec, lig, simPar, table, &h_trafoLig, &h_potLig);

			// rotate forces back to global frame here, if it might become neccessary.

//			//DEBUG
//			if(count++ == 18494 && true) {
//				h_potLigNL.resetHostData();
//				IP.h_NLPotForce(grid, rec, lig, simPar, table, &h_trafoLig, &h_potLigNL);
//				h_potLigNL.printEl(lig->numAtoms(), 0);
////				return;
//			}

			TF.h_RedEnForce(lig, dof, rotMat, &h_potLig, NULL, engrad);

			nvtxRangePop();

		}

		item->setReady();
	}

#else
	using namespace std;
	int precisionSetting = cerr.precision( );
	ios::fmtflags flagSettings = cerr.flags();
	cerr.setf(ios::scientific | ios::showpos);
	cerr.precision(5);
	/* Initialize Input/Output Buffers */

	//DEBUG
//	Comp5_HD<float, HOSTONLY> h_potLigNL(_atomBufferSize);
//	h_potLigNL.initHost();

	/* Create local Transformer and Interpolator object */
	Transformer TF;
	Interpolator IP;

	/* Create SimParam obj */
	/* ToDo: object of this type should be added to the dataMngt in global and per-client manner */
	SimParam simPar;
	simPar.dielec = variable;				/** type of dielectric constant */
	simPar.epsilon = 15;					/** dielectric constant */
	simPar.ffelec = FELEC/simPar.epsilon;	/** precomputed factor felec/epsilon */
//	simPar.useSwi = false;					/** using switching potential */
//	simPar.swiOn = 0;						/** min. switching potential distance */
//	simPar.swiOff= 0;						/** max. switching potential distance */
	simPar.useRecGrad = false;				/** using Receptor gradients */
	simPar.usePot = true;					/** use Potential grid */

	/* Only one item available */
	/* get a workItem from the queue */
	WorkerItem* item;
	nvtxRangePushA("CPUWorker BlockedRemove");
	item = _queue.removeCondBlocked();
	nvtxRangePop();

	assert(item->size() > 0);

	/* each DOF is processed sequentially */

	const Protein* rec = getProtein(item->globRecId());
	const Protein* lig = getProtein(item->globLigId());
	const GridUnion* grid = getGridUnion(item->globGridId());
	const AttrParamTable* table = getParamTable();

	unsigned atomBufferSize = lig->numAtoms();
	Comp3_HD<float, HOSTONLY> h_trafoLig(atomBufferSize);
	h_trafoLig.initHost();

	Comp3_HD<float, HOST_PINNED> d_trafoLig(atomBufferSize);
	d_trafoLig.initDevice();
	d_trafoLig.set_h_x(h_trafoLig.h_x());
	d_trafoLig.set_h_y(h_trafoLig.h_y());
	d_trafoLig.set_h_z(h_trafoLig.h_z());

	Comp5_HD<float, HOSTONLY> h_potLig(atomBufferSize);
	h_potLig.initHost();

	Comp5_HD<float, HOST_PINNED> d_potLig(atomBufferSize);
	d_potLig.initDevice();
	d_potLig.set_h_x(h_potLig.h_x());
	d_potLig.set_h_y(h_potLig.h_y());
	d_potLig.set_h_z(h_potLig.h_z());
	d_potLig.set_h_v(h_potLig.h_v());
	d_potLig.set_h_w(h_potLig.h_w());

	Timer timer;
	cudaEvent_t start, stop;
	cudaVerify(cudaEventCreate(&start));
	cudaVerify(cudaEventCreate(&stop));

	double tDOF2Coord, tPotForce, tNL, tForce2Grad ;
	tDOF2Coord = tPotForce = tNL = tForce2Grad = 0.0;

	for (int i = 0; i < item->size(); ++i) {
		const DOF& dof = item->DOFBuffer()[i];
		EnGrad& engrad = item->EnGradBuffer()[i];
		RotMatf rotMat;
		initTimer(&timer);
		core::euler2rotmat(dof.ang.x, dof.ang.y, dof.ang.z, rotMat);

		initTimer(&timer);
		TF.h_DOF2Pos(lig, dof, rotMat, &h_trafoLig);
		tDOF2Coord += getTimer(&timer);

		initTimer(&timer);
		IP.h_PotForce(grid,lig, &h_trafoLig, &h_potLig);
		tPotForce += getTimer(&timer);

		initTimer(&timer);
		IP.h_NLPotForce(grid, rec, lig, simPar, table, &h_trafoLig, &h_potLig);
		tNL += getTimer(&timer);

		initTimer(&timer);
		TF.h_RedEnForce(lig, dof, rotMat, &h_potLig, NULL, engrad);
		tForce2Grad += getTimer(&timer);

	}

	std::cerr	<< "TIME" << "\t"
				<< tDOF2Coord << "\t"
				<< tPotForce << "\t"
				<< tNL << "\t"
				<< tForce2Grad << std::endl;

	cerr.precision(precisionSetting);
	cerr.flags(flagSettings);

	item->setReady();
	/* makes the worker sleeping */
	item = _queue.removeCondBlocked();

#endif
	global_log->info() << std::setw(5) << " " << "CPUWorker " << _id <<" terminated" << std::endl;


}


/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/
