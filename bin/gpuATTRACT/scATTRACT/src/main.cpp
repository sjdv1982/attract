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

#include <memory>
#include <string>
#include <fstream>
#include <nvToolsExt.h>

#include <AttractServer>
#include "asUtils/Logger.h"
#include "asClient/DOFTransform.h"
#include "asClient/cudaArchCheck.h"
#include "asUtils/timer.h"

#include "tclap/CmdLine.h"

using namespace std;

namespace score {
Log::Logger* _log;
}


void init_logger( bool use_file = true) {
	using namespace Log;
#ifndef NDEBUG
	logLevel level = Debug;
#else
	logLevel level = Info;
#endif
	if (use_file) {
		string filename = "ScoringClient.log";
		score::_log = new Logger(level, filename.substr(0,filename.size()-4));
	} else {
		score::_log = new Logger(level, &(std::cerr));
	}
}

void printResultsScore(unsigned numDOFs, as::DOF* dofs, as::EnGrad* enGrads);

int main (int argc, char *argv[]) {

	/* initialize Logger */
	bool use_file = true;
	init_logger(use_file);
	unique_ptr<Log::Logger> log(score::_log);

	/* required variables */
	string dofName;

	/* optional variables */
	string gridName;
	string ligName;
	string recName;
	string paramsName;
	string recGridAlphabetName;

	int numCPUs;
	vector<int> devices;
	int chunkSize;
	int maxItemsPerSubmit;

	int numToConsider;
	int whichToTrack;
	/* catch command line exceptions */
	try {

		/* print argv */
		log->info() << "Client starts with command: ";
		std::vector<std::string> arguments(argv , argv + argc);
		for (auto string : arguments) { *log << string << " ";}
		*log << endl;

		/* description of the application */
		TCLAP::CmdLine cmd("An ATTRACT client that performs scoring.", ' ', "1.1");

		/* define required arguments */
		TCLAP::ValueArg<string> dofArg("","dof","",true,"Structure (DOF) file","*.dat", cmd);

		/* define optional arguments */
		TCLAP::ValueArg<string> recArg("r","receptor-pdb","pdb-file name of receptor. (Default: receptorr.pdb)", false,"receptorr.pdb","*.pdb", cmd);
		TCLAP::ValueArg<string> ligArg("l","ligand-pdb","pdb-file name of ligand. (Default: ligandr.pdb)", false, "ligandr.pdb","*.pdb", cmd);
		TCLAP::ValueArg<string> gridArg("g","grid","Receptor grid file. (Default: receptorgrid.grid)",false, "receptorgrid.grid","*.grid", cmd);
		TCLAP::ValueArg<string> paramArg("p","par","Attract parameter file. (Default: attract.par)",false,"attract.par","*.par", cmd);
		TCLAP::ValueArg<string> gridAlphabet("a","receptor-alphabet","Receptor grid alphabet file.",false,"","*.alphabet", cmd);

		TCLAP::ValueArg<int> cpusArg("c","cpus","Number of CPU threads to be used. (Default: 0)", false, 0, "int", cmd);

		int numDevicesAvailable = 0; CUDA_CHECK(cudaGetDeviceCount(&numDevicesAvailable));
		vector<int> allowedValues(numDevicesAvailable); iota(allowedValues.begin(), allowedValues.end(), 0);
		TCLAP::ValuesConstraint<int> vc(allowedValues);
		TCLAP::MultiArg<int> deviceArg("d","device","Device ID of GPU to be used.", false, &vc, cmd);
		TCLAP::ValueArg<int> chunkSizeArg("","chunkSize", "Number of concurrently processed structures", false, 1000, "int", cmd);
		TCLAP::ValueArg<int> maxItemsArg("","maxItems", "Max. number of item per submit", false, 10000, "int", cmd);

		TCLAP::ValueArg<int> num2ConsiderArg("","num", "Number of configurations to consider (1 - num). (Default: All)", false, -1, "int", cmd);
		TCLAP::ValueArg<int> which2TrackArg("","focusOn", "Condider only this configuration. (Default: -1)", false, -1, "int", cmd);


		// parse cmd-line input
		cmd.parse(argc, argv);

		/* Assigne parsed values */
		recName 	= recArg.getValue();
		ligName 	= ligArg.getValue();
		gridName 	= gridArg.getValue();
		paramsName 	= paramArg.getValue();
		dofName 	= dofArg.getValue();
		devices 	= deviceArg.getValue();
		numCPUs 	= cpusArg.getValue();
		chunkSize 	= chunkSizeArg.getValue();
		maxItemsPerSubmit = maxItemsArg.getValue();
		numToConsider = num2ConsiderArg.getValue();
		whichToTrack = which2TrackArg.getValue();
		recGridAlphabetName = gridAlphabet.getValue();

	} catch (TCLAP::ArgException &e){
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
	}


	log->info() << "recName=" << recName 		<< endl;
	log->info() << "ligName=" << ligName 		<< endl;
	log->info() << "gridName=" << gridName 		<< endl;
	log->info() << "parName=" << paramsName 	<< endl;
	log->info() << "dofName=" << dofName	 	<< endl;
	log->info() << "recGridAlphabetName" << recGridAlphabetName << endl;
	log->info() << "numCPUs=" << numCPUs 		<< endl;
	log->info() << "devices=[ "; for (auto device : devices) *log << device << " "; *log << "]"<<  endl;
	log->info() << "chunkSize=" << chunkSize 	<< endl;
	log->info() << "maxItems=" << maxItemsPerSubmit	 <<  endl;
	log->info() << "numToConsider=" << numToConsider << endl;
	log->info() << "whichToTrack=" << whichToTrack << endl;

	if (devices.size() > 0) {
		/* Check Compute Capability of devices */
		try {
			asClient::checkComputeCapability();
		} catch (std::exception& e) {
			cerr << "Error: " << e.what() << endl;
			exit(EXIT_FAILURE);
		}
	}

	/* read dof header */
	std::vector<asUtils::Vec3f> pivots;
	bool autoPivot;
	bool centered_receptor, centered_ligands;
	asDB::readDOFHeader(dofName, pivots, autoPivot, centered_receptor, centered_ligands);

	/* check file. only a receptor-ligand pair (no multi-bodies!) is allowed */
	if(!autoPivot && pivots.size() > 2) {
		log->error() << "DOF-file contains defintions for more than two molecules. Multi-body docking is not supported." << endl;
		exit(EXIT_FAILURE);
	}

	/* read dofs */
	std::vector<std::vector<as::DOF>> DOF_molecules;
	asDB::readDOFFromFile(dofName, DOF_molecules);
	/* check file. only one receptor-ligand pair (no multi-bodies!) is allowed */
	if(DOF_molecules.size() != 2) {
		log->error() << "DOF-file contains defintions for more than two molecules. Multi-body docking is not supported." << endl;
		exit(EXIT_FAILURE);
	}

	/* shrink number of dofs artificially */
	if (numToConsider >= 0 || whichToTrack >= 0) {
		if (whichToTrack >= 0) {
			DOF_molecules[1][0] = DOF_molecules[1][whichToTrack];
			DOF_molecules[1].resize(1);
			DOF_molecules[0].resize(1);
		} else {
			DOF_molecules[1].resize(numToConsider);
			DOF_molecules[0].resize(numToConsider);
		}
	}

	/* initialize the scoring server: mngt(numItems, chunkSize, deviceBufferSize )*/
	/* only a maximum number of items per request is allowed since too many items introduce a large overhead.
	 * To circumvent this problem, we split into multiple requests */

	int ligandSize = asDB::readProteinSizeFromPDB(ligName);
	int deviceBufferSize = ligandSize*chunkSize;
	int numDofs = DOF_molecules[0].size();
	int numItems = (numDofs + chunkSize - 1) /chunkSize;

	log->info() << "deviceBufferSize=" << deviceBufferSize << endl;
	log->info() << "ligandSize=" << ligandSize << endl;
	log->info() << "numDofs=" << numDofs << endl;
	log->info() << "numItems=" << numItems << endl;

	int numSubmitIter = 1; // number of request submission
	int strucPerSubmit = numDofs;
	int strucPerSubmitLast = 0;
	int wait_ms = 2*180000; // max wait for pulling request in ms

	if (numItems > maxItemsPerSubmit) {
		/* choose items per submit so that it is multiple of the chunkSize */
		log->info() << "DOF array is splitted in multiple requests" << endl;

		strucPerSubmit = maxItemsPerSubmit*chunkSize;
		numSubmitIter = (numDofs + strucPerSubmit - 1) / strucPerSubmit;
		strucPerSubmitLast = numDofs - (numSubmitIter - 1)*strucPerSubmit;

		log->info() << "strucPerSubmit=" << strucPerSubmit << "numSubmitIter=" << numSubmitIter
				<< "strucPerSubmitLast=" << strucPerSubmitLast << endl;
		numItems = 2*maxItemsPerSubmit;

	}

	as::ServerManagement server(numItems, chunkSize, deviceBufferSize);


	/* load proteins and grid, and get a handle to it*/
	const int clientId = 0; // by specifing a client id, we may remove all added data by the by a call to removeClient(id)
	int ligId, recId, gridId;
	recId = server.addProtein(clientId, recName);
	ligId = server.addProtein(clientId ,ligName);
	gridId = server.addGridUnion(clientId, gridName);
	server.addParamTable(paramsName);

//	server.getProtein(recId)->print();

	if(autoPivot) {
		if (!pivots.empty()) {
			log->error() << "Auto pivot specified, but explicitly definined pivots available. (File "<< dofName << ")" << endl;
			exit(EXIT_FAILURE);
		}
		pivots.push_back(server.getProtein(recId)->pivot());
		pivots.push_back(server.getProtein(ligId)->pivot());
	} else {
		if (pivots.size() != 2) {
			log->error() << "No auto pivot specified, but number of definined pivots is incorrect. (File "<< dofName << ")" << endl;
			exit(EXIT_FAILURE);
		}
		server.getProtein(recId)->pivotize(pivots[0]);
		server.getProtein(ligId)->pivotize(pivots[1]);
	}

	log->info() << "pivots= "; for (auto pivot : pivots) *log << pivot << ", "; *log << endl;

	/* apply receptor grid mapping for ligand */
	if (!recGridAlphabetName.empty()) {
		std::vector<unsigned> mapVec = asDB::readGridAlphabetFromFile(recGridAlphabetName);
		as::TypeMap typeMap = as::createTypeMapFromVector(mapVec);
		as::Protein* prot = server.getProtein(ligId);
		as::applyDefaultMapping(prot->numAtoms(), prot->type(), prot->type());
		as::applyMapping(typeMap, prot->numAtoms(), prot->type(), prot->mappedTypes());
	} else {
		log->warning() << "No grid alphabet specified. Applying default mapping." << endl;
		as::Protein* prot = server.getProtein(ligId);
		as::applyDefaultMapping(prot->numAtoms(), prot->type(), prot->type());
		as::applyDefaultMapping(prot->numAtoms(), prot->type(), prot->mappedTypes());
	}

	if (true){
		as::GridUnion* grid = server.getGridUnion(gridId);
		asUtils::Vec3f& pivot = pivots[0];
		float3 pos;

		pos = grid->innerGrid()->pos();
		pos.x -= pivot[0]; pos.y -= pivot[1]; pos.z -= pivot[2];
		grid->innerGrid()->setPos(pos);
//		as::print(grid->innerGrid(), 0, 0 , 0,0, 25, 30, 25, 30, 25, 30);

		pos = grid->outerGrid()->pos();
		pos.x -= pivot[0]; pos.y -= pivot[1]; pos.z -= pivot[2];
		grid->outerGrid()->setPos(pos);

		pos = grid->NLgrid()->pos();
		pos.x -= pivot[0]; pos.y -= pivot[1]; pos.z -= pivot[2];
		grid->NLgrid()->setPos(pos);
	}

//	server.getProtein(recId)->print();

	/* initialize CPU workers */
	for(int i = 0; i < numCPUs; ++i) {
		server.addCPUWorker();
	}


	/* transform ligand dofs assuming that the receptor is always centered in the origin */
	asClient::transformDOF_glob2rec(DOF_molecules[0], DOF_molecules[1], pivots[0], pivots[1], centered_receptor, centered_ligands);
	as::DOF* dofBuffer = DOF_molecules[1].data();

	/* Allocate result buffer */
	vector<as::EnGrad> GPU_enGrad(numDofs);
	vector<as::EnGrad> CPU_enGrad(numDofs);

	as::EnGrad* GPU_enGradBuffer = GPU_enGrad.data();
    as::EnGrad* CPU_enGradBuffer = CPU_enGrad.data();


	/* Perform calculations */
	asUtils::Timer timer;


	//DEBUG
//	AServ::core::nlCount = 0;

	/* CPUWorker */
	if (numCPUs > 0) {

		int reqIds[2] = {-1, -1};

		asUtils::initTimer(&timer);

		int count = 0;

		nvtxRangePushA("Main Submit");
		while (reqIds[0] == -1 && count < wait_ms) {
			reqIds[0] = server.submitRequest(dofBuffer, strucPerSubmit, gridId, recId, ligId, as::Request::CPU);
			if (reqIds[0] >= 0) {
				break;

			}
			++count;
			std::this_thread::sleep_for(std::chrono::milliseconds(1));

		}
		nvtxRangePop();
		if (count >= 1) {
			log->error() << "Error: Submitting request." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		log->info() << "Submit " << 0 << " Id " << reqIds[0] <<  endl;

		std::swap(reqIds[0], reqIds[1]);
		/* submit next requests */
		for (int i = 1; i < numSubmitIter; ++i) {
			nvtxRangePushA("Main Submit");
			while (reqIds[0] == -1 && count < wait_ms) {
				reqIds[0] = server.submitRequest(dofBuffer + i*strucPerSubmit, i == (numSubmitIter -1) ?  strucPerSubmitLast : strucPerSubmit, gridId, recId, ligId, as::Request::CPU);
				if (reqIds[0] >= 0) {
					break;

				}
				++count;
				std::this_thread::sleep_for(std::chrono::milliseconds(1));

			}
			nvtxRangePop();
			if (count >= 1) {
				log->error() << "Error: Submitting request." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			log->info() << "Submit " << i << " Id " << reqIds[0] <<  endl;

			while ((server.pullRequest(reqIds[1], CPU_enGradBuffer + (i-1)*strucPerSubmit) != as::Dispatcher::ready) && count++ < wait_ms) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
			if (count >= wait_ms) {
				log->error() << "Error: Pulling request." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			log->info() << "Pull " << i-1 << " Id " << reqIds[1] <<  endl;

			reqIds[1] = -1;
			count = 0;

			std::swap(reqIds[0], reqIds[1]);

		}

		/* Pull for last request */
		while ((server.pullRequest(reqIds[1], CPU_enGradBuffer + (numSubmitIter-1)*strucPerSubmit) != as::Dispatcher::ready) && count++ < wait_ms) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}

		asUtils::getTimerAndPrint(&timer, "CPU");

		if (count >= wait_ms) {
			log->error() << "Error: Pulling request." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		log->info() << "Pull " << numSubmitIter-1 << " Id " << reqIds[1] <<  endl;


	}


	/* GPUWorker */
	if (devices.size() > 0) {

		/* initialize GPU workers and add/transfere the necessary data to the devices.
		 * Note, in general, not all data need to be available at each GPU */
		for (auto deviceId : devices) {
			server.addGPUWorker(deviceId);
			server.attachProteinToDevice(recId, deviceId);
			server.attachProteinToDevice(ligId, deviceId);
			server.attachGridUnionToDevice(gridId, deviceId);
			server.attachParamTableToDevice(deviceId);
		}

		/* after the data transfere, we need to call updateDeviceIDLookup() to let the dispatcher know
		 * which data is on which device */
		if (devices.size() > 0) {
			server.updateDeviceIDLookup();
		}

		int reqIds[2] = {-1, -1};

		asUtils::initTimer(&timer);

		int count = 0;

		nvtxRangePushA("Main Submit");
		while (reqIds[0] == -1 && count < wait_ms) {
			reqIds[0] = server.submitRequest(dofBuffer, strucPerSubmit, gridId, recId, ligId, as::Request::GPU);
			if (reqIds[0] >= 0) {
				break;

			}
			++count;
			std::this_thread::sleep_for(std::chrono::milliseconds(1));

		}
		nvtxRangePop();
		if (count >= 1) {
			log->error() << "Error: Submitting request." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		log->info() << "Submit " << 0 << " Id " << reqIds[0] <<  endl;

		std::swap(reqIds[0], reqIds[1]);
		/* submit next requests */
		for (int i = 1; i < numSubmitIter; ++i) {
			nvtxRangePushA("Main Submit");
			while (reqIds[0] == -1 && count < wait_ms) {
				reqIds[0] = server.submitRequest(dofBuffer + i*strucPerSubmit, i == (numSubmitIter -1) ?  strucPerSubmitLast : strucPerSubmit, gridId, recId, ligId, as::Request::GPU);
				if (reqIds[0] >= 0) {
					break;

				}
				++count;
				std::this_thread::sleep_for(std::chrono::milliseconds(1));

			}
			nvtxRangePop();
			if (count >= 1) {
				log->error() << "Error: Submitting request." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			log->info() << "Submit " << i << " Id " << reqIds[0] <<  endl;

			while ((server.pullRequest(reqIds[1], GPU_enGradBuffer + (i-1)*strucPerSubmit) != as::Dispatcher::ready) && count++ < wait_ms) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
			if (count >= wait_ms) {
				log->error() << "Error: Pulling request." << std::endl;
				std::exit(EXIT_FAILURE);
			}

			log->info() << "Pull " << i-1 << " Id " << reqIds[1] <<  endl;

			reqIds[1] = -1;
			count = 0;

			std::swap(reqIds[0], reqIds[1]);

		}

		/* Pull for last request */
		while ((server.pullRequest(reqIds[1], GPU_enGradBuffer + (numSubmitIter-1)*strucPerSubmit) != as::Dispatcher::ready) && count++ < wait_ms) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}

		asUtils::getTimerAndPrint(&timer, "GPU");

		if (count >= wait_ms) {
			log->error() << "Error: Pulling request." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		log->info() << "Pull " << numSubmitIter-1 << " Id " << reqIds[1] <<  endl;

	}


	if (numCPUs > 0 != devices.size() > 0) {
		if (numCPUs > 0) {
			printResultsScore(numDofs, dofBuffer, CPU_enGradBuffer);
		}
		if (devices.size() > 0) {
			printResultsScore(numDofs, dofBuffer, GPU_enGradBuffer);
		}
	}



	return 0;
}

void printResultsScore(unsigned numDOFs, as::DOF* dofs, as::EnGrad* enGrads)
{
	using namespace std;

	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();

	for (int i = 0; i < numDOFs; ++i) {
		const as::EnGrad& enGrad = enGrads[i];
		const as::DOF& dof = dofs[i];
		cout.setf(ios::scientific);
		cout.precision(8);
		cout << " Energy: " << enGrad.E_VdW + enGrad.E_El << endl;
		cout.unsetf(ios::scientific);

		cout.setf(ios::fixed);
		cout.precision(3);
		cout << setw(12) << enGrad.E_VdW << setw(12) << enGrad.E_El << endl;
		cout.unsetf(ios::fixed);

		cout.setf(ios::scientific);
		cout.precision(8);
		int width = 20;
		cout << " Gradients: "
				<< setw(width) << enGrad.ang.x  << setw(width) << enGrad.ang.y  << setw(width) << enGrad.ang.z
				<< setw(width) << enGrad.pos.x  << setw(width) << enGrad.pos.y  << setw(width) << enGrad.pos.z  << endl;
		cout.unsetf(ios::scientific);

	}

	cout.precision(precisionSetting);
	cout.flags(flagSettings);
}
