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

#include <iostream>
#include <sstream>
#include <list>
#include <cassert>
#include <random>
#include <cmath>
#include <memory>
#include <string>
#include <nvToolsExt.h>
#include <cuda_runtime.h>

#include <AttractServer>
#include "asUtils/Logger.h"
#include "asClient/DOFTransform.h"
#include "asClient/cudaArchCheck.h"
#include "asUtils/timer.h"

#include "tclap/CmdLine.h"

using namespace std;

namespace mca {
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
		string filename = "mcATTRACT.log";
		mca::_log = new Logger(level, filename.substr(0,filename.size()-4));
	} else {
		mca::_log = new Logger(level, &(std::cerr));
	}
}

/*
 ** @brief: creates a random rotation axis that is uniformly distributed on a spheres surface
 *
 ** @input: rand0 & rand1 are random numbers in the range [0, 1].
 *
 ** @output: axis[3], vector of length 1.0 .
 */
double PIx2 = 2*M_PI;
inline void randAxis(const double& rand0, const double& rand1, double axis[3]) {
	double phi = PIx2 * rand1;
	double theta = std::acos(2*rand0 - 1);
	double stheta = std::sin(theta);

	axis[0] = stheta*cos(phi);
	axis[1] = stheta*sin(phi);
	axis[2] = cos(theta);
}

/*
 ** @brief: creates a rotation matrix given an axis and an angle
 */
inline void rotMatFromAxisAng(const double axis[3], const double ang, asUtils::RotMatd& rotMat) {
	double c = std::cos(ang);
    double s = std::sin(ang);
    double t = 1.0 - c;

    rotMat[0] = c + axis[0]*axis[0]*t;
    rotMat[4] = c + axis[1]*axis[1]*t;
    rotMat[8] = c + axis[2]*axis[2]*t;


    double tmp1 = axis[0]*axis[1]*t;
    double tmp2 = axis[2]*s;
    rotMat[1] = tmp1 - tmp2;
    rotMat[3] = tmp1 + tmp2;

    tmp1 = axis[0]*axis[2]*t;
    tmp2 = axis[1]*s;
    rotMat[2] = tmp1 + tmp2;
    rotMat[6] = tmp1 - tmp2;


    tmp1 = axis[1]*axis[2]*t;
    tmp2 = axis[0]*s;
    rotMat[5] = tmp1 - tmp2;
    rotMat[7] = tmp1 + tmp2;
}


static double maxDist;
static double maxAng;
static double kT;

/* Initialize random number generators */
static std::default_random_engine generator;
static std::uniform_real_distribution<double> distribution(0.0, 1.0);

inline void randomStep (const as::DOF& oldDOF, as::DOF& newDOF)
{
//	int numdis = 20;
//	static int count = 0;
//	if (++count < numdis) {
//		std::cout << "#" << count << std::endl;
//	}
	/* get 6 random numbers: 3 pos. + 3 ang. */
	double r[6];
	for (unsigned k = 0; k < 6; ++k) {
		r[k] = distribution(generator);
	}

//	if (count < numdis) {
//		std::cout << "rand "<< r[0] << " "<< r[1] << " "<< r[2] << " "<< r[3] << " "<< r[4] << " "<< r[5] << std::endl;
//	}

	/********** Apply random rotation **********/

	/* create random axis */
	double rAxis[3];
	randAxis(r[0], r[1], rAxis);
//	if (count < numdis) {
//		std::cout << "rAxis "<< r[0] << " "<< r[1] << " "<< r[2] << " "<< r[3] << " "<< r[4] << " "<< r[5] << std::endl;
//	}

	/* calculate rotation matrix according to axis and angle */
	double ang = (2.0*r[2]-1.0)*maxAng; // rad !!!
	asUtils::RotMatd randRotMat;

	rotMatFromAxisAng(rAxis, ang, randRotMat);
//	if (count < numdis) {
//		std::cout << "randRotMat " << randRotMat << endl;
//	}

	/* calculate current rotation matrix accoring to euler angles */
	asUtils::RotMatd currRotMat;
	double phi, ssi, rot;
	phi = (double) oldDOF.ang.x;
	ssi = (double) oldDOF.ang.y;
	rot = (double) oldDOF.ang.z;

	asClient::euler2rotmat(phi, ssi, rot, currRotMat);
//	if (count < numdis) {
//		std::cout << "currRotMat " << currRotMat << endl;
//	}



	/* Multiply matrices and retain new angels */
	asUtils::RotMatd newRotmat;
	newRotmat = randRotMat*currRotMat;
	asClient::rotmat2euler(newRotmat, phi, ssi, rot);
	newDOF.ang.x = (float)phi;
	newDOF.ang.y = (float)ssi;
	newDOF.ang.z = (float)rot;


	/********** Apply random displacement **********/

	double dx = (2.0*r[3] - 1.0)*maxDist;
	double dy = (2.0*r[4] - 1.0)*maxDist;
	double dz = (2.0*r[5] - 1.0)*maxDist;

//	cout << "#" << count << endl;
//	cout << rAxis[0] << " " << rAxis[1] << " " << rAxis[2] << " " << ang * 180.0 / M_PI << endl;
//	cout << dx << " " << dy << " " << dz << endl;
//	count++;

	newDOF.pos.x = oldDOF.pos.x + dx;
	newDOF.pos.y = oldDOF.pos.y + dy;
	newDOF.pos.z = oldDOF.pos.z + dz;

//	if (count < numdis) {
//		std::cout << "old " << oldDOF << std::endl;
//		std::cout << "new " << newDOF << std::endl;
//	}
}


void MC_accept(as::DOF& oldDOF, as::EnGrad& oldEnGrad, as::DOF &newDOF, as::EnGrad& newEnGrad) {
	float newEnergy = newEnGrad.E_El + newEnGrad.E_VdW;
	float oldEnergy = oldEnGrad.E_El + oldEnGrad.E_VdW;

//	bool accepted = false;
	/* Metropolis Criterion */
	if (newEnergy <= oldEnergy) {
		oldEnGrad = newEnGrad;
		oldDOF = newDOF;
//		accepted = true;
	} else {
		double r = distribution(generator);
		if (r < std::exp(-(newEnergy - oldEnergy)/kT)) {
			oldEnGrad = newEnGrad;
			oldDOF = newDOF;
//			accepted = true;
		}
	}

//	int numdis = 20;
//	static int count = 0;
//	if (++count < numdis) {
//		std::cout << "accepted=" << accepted << " newE="<< newEnergy << " oldE=" << oldEnergy << std::endl;
//	}
}

unsigned readProteinSizeFromFile (std::string name) {
	using namespace std;

	string filename = name;

	ifstream file(filename.c_str(), ios::in | ios::binary);
	ProteinDesc desc;

	if (!file.read((char*)&desc.numAtoms, sizeof(unsigned))) {
		cerr << "Error read: numAtoms" << endl;
		exit(EXIT_FAILURE);
	}
	return desc.numAtoms;
}

/* printing results to stdout */
void printResultsOutput(unsigned numDofs, as::DOF* dofs, as::EnGrad* enGrads, std::vector<asUtils::Vec3f>& pivots)
{
	using namespace std;

	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();
	cout.setf(ios::showpoint);
	cout.precision(6);

	asUtils::Vec3f pivot_diff = pivots[0] - pivots[1];

	/* print header */
	cout << "#pivot 1 " << pivots[0][0] << " " << pivots[0][1] << " " << pivots[0][2] << " " << endl;
	cout << "#pivot 2 " << pivots[1][0] << " " << pivots[1][1] << " " << pivots[1][2] << " " << endl;
	cout << "#centered receptor: false" << endl;
	cout << "#centered ligands: false" << endl;
	for (unsigned i = 0; i < numDofs; ++i) {
		const as::EnGrad& enGrad = enGrads[i];
		const as::DOF& dof = dofs[i];
		cout << "#"<< i+1 << endl;
		cout << "## Energy: " << enGrad.E_VdW + enGrad.E_El << endl;
		cout << "## " << enGrad.E_VdW << " " << enGrad.E_El << endl;
		cout << 0.0 << " " << 0.0 << " " << 0.0 << " "
			 << 0.0 << " " << 0.0 << " " << 0.0 << endl;
		cout << dof.ang.x << " " << dof.ang.y << " " << dof.ang.z << " "
			 << dof.pos.x + pivot_diff[0]<< " " << dof.pos.y + pivot_diff[1] << " " << dof.pos.z + pivot_diff[2] << endl;
	}

	cout.precision(precisionSetting);
	cout.flags(flagSettings);
}


int main (int argc, char *argv[]) {
	using namespace std;

	/* initialize Logger */
	bool use_file = true;
	init_logger(use_file);
	unique_ptr<Log::Logger> log(mca::_log);

	/* required variables */
	string dofName;

	/* optional variables */
	string gridName;
	string ligName;
	string recName;
	string paramsName;
	string recGridAlphabetName;
//	static double maxDist;
//	static double maxAng;	They are set as well;
//	static double kT;
	unsigned numCPUs;
	unsigned numIter;
	unsigned chunkSize;
	vector<int> devices;

	/* catch command line exceptions */
	try {

		/* print argv */
		log->info() << "Client starts with command: ";
		std::vector<std::string> arguments(argv , argv + argc);
		for (auto string : arguments) { *log << string << " ";}
		*log << endl;

		/* description of the application */
		TCLAP::CmdLine cmd("An ATTRACT client that performs energy minimization by a Monte Carlo search.", ' ', "1.1");

		/* define required arguments */
		TCLAP::ValueArg<string> dofArg("","dof","",true,"Structure (DOF) file","*.dat", cmd);

		/* define optional arguments */
		TCLAP::ValueArg<string> recArg("r","receptor-pdb","pdb-file name of receptor. (Default: receptorr.pdb)", false,"receptorr.pdb","*.pdb", cmd);
		TCLAP::ValueArg<string> ligArg("l","ligand-pdb","pdb-file name of ligand. (Default: ligandr.pdb)", false, "ligandr.pdb","*.pdb", cmd);
		TCLAP::ValueArg<string> gridArg("g","grid","Receptor grid file. (Default: receptorgrid.grid)",false, "receptorgrid.grid","*.grid", cmd);
		TCLAP::ValueArg<string> paramArg("p","par","Attract parameter file. (Default: attract.par)",false,"attract.par","*.par", cmd);
		TCLAP::ValueArg<string> gridAlphabet("a","receptor-alphabet","Receptor grid alphabet file.",false,"","*.alphabet", cmd);

		TCLAP::ValueArg<unsigned> cpusArg("c","cpus","Number of CPU threads to be used. (Default: 0)", false, 0, "uint");
		TCLAP::ValueArg<unsigned> numIterArg("","iter","Number Monte Carlo iterations. (Default: 50)", false, 50, "uint", cmd);
		TCLAP::ValueArg<double> maxDistArg("","maxDist","Maximum translational displacement (A). (Default: 1.0A)", false, 1.0, "int", cmd);
		TCLAP::ValueArg<double> maxAngArg("","maxAng","Maximum rotational displacement (deg). (Default: 3.0deg)", false, 3.0, "int", cmd);
		TCLAP::ValueArg<double> kTArg("","kT","Monte Carlo temperature. (Default: 10.0)", false, 10.0, "int", cmd);
		int numDevicesAvailable = 0; CUDA_CHECK(cudaGetDeviceCount(&numDevicesAvailable));;
		vector<int> allowedValues(numDevicesAvailable); iota(allowedValues.begin(), allowedValues.end(), 0);
		TCLAP::ValuesConstraint<int> vc(allowedValues);
		TCLAP::MultiArg<int> deviceArg("d","device","Device ID of serverMode to be used.", false, &vc);
		TCLAP::ValueArg<unsigned> chunkSizeArg("","chunkSize", "Number of concurrently processed structures", false, 5000, "uint", cmd);

		cmd.xorAdd(cpusArg, deviceArg);

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
		numIter		= numIterArg.getValue();
		maxDist		= maxDistArg.getValue();
		maxAng		= maxAngArg.getValue();
		kT 			= kTArg.getValue();
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
	log->info() << "numIter=" << numIter 		<< endl;
	log->info() << "maxDist=" << maxDist		<< endl;
	log->info() << "maxAng="  << maxAng			<< endl;
	log->info() << "kT="  	  << kT				<< endl;
	log->info() << "chunkSize=" << chunkSize 	<< endl;

	if (devices.size() > 0) {
		/* Check Compute Capability of devices */
		try {
			asClient::checkComputeCapability();
		} catch (std::exception& e) {
			cerr << "Error: " << e.what() << endl;
			exit(EXIT_FAILURE);
		}
	}

	/* convert degrees to rad */
	maxAng = maxAng * M_PI / 180.0;

	/* check if cpu or gpu is used */
	as::Request::useMode_t serverMode = as::Request::unspecified;
	if(numCPUs > 0) {
		serverMode = as::Request::CPU;
	} else if (devices.size() > 0) {
		serverMode = as::Request::GPU;
	} else {
		log->error() << "Neither CPU nor GPU is specified. This state should not happen." << endl;
		exit(EXIT_FAILURE);
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

	/*
	 * initialize the scoring server: mngt(numItems, chunkSize, deviceBufferSize )
	 * only a maximum number of items per request is allowed since too many items introduce a large overhead
	 */

	unsigned ligandSize = asDB::readProteinSizeFromPDB(ligName);
	unsigned deviceBufferSize = ligandSize*chunkSize;
	unsigned numDofs = DOF_molecules[0].size();
	const unsigned numItems = 4*(((unsigned)ceil((double)numDofs/2.0) + chunkSize - 1) / chunkSize);

	log->info() << "ligandSize=" 		<< ligandSize 		<< endl;
	log->info() << "deviceBufferSize=" 	<< deviceBufferSize << endl;
	log->info() << "numDofs=" 			<< numDofs 			<< endl;
	log->info() << "numItems=" 			<< numItems 		<< endl;

	/* check if there are too two many items per request */
	const int maxItemsPerSubmit = 10000;
	if (numItems/4 > maxItemsPerSubmit) {
		log->error() << "Too many items per request. Increase chunkSize" << endl;
		exit(EXIT_FAILURE);
	}

	as::ServerManagement server(numItems, chunkSize, deviceBufferSize);

	/* load proteins and grid, and get a handle to it*/
	const int clientId = 0; // by specifing a client id, we may remove all added data by the by a call to removeClient(id)
	int ligId, recId, gridId;
	recId = server.addProtein(clientId, recName);
	ligId = server.addProtein(clientId ,ligName);
	gridId = server.addGridUnion(clientId, gridName);
	server.addParamTable(paramsName);

	/* parse or get pivots. only two pivots/molecules are allowed */
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

	/* transform ligand dofs assuming that the receptor is always centered in the origin */
	asClient::transformDOF_glob2rec(DOF_molecules[0], DOF_molecules[1], pivots[0], pivots[1], centered_receptor, centered_ligands);

	/* adapt grid locations according to receptor pivot (actually, I don't get why we need this) */
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

	/* initialize CPU workers if any*/
	for(unsigned i = 0; i < numCPUs; ++i) {
		server.addCPUWorker();
	}

	/* initialize devices and GPU workers if any */
	for (unsigned i = 0; i < devices.size(); ++i) {
		server.attachProteinToDevice(recId, devices[i]);
		server.attachProteinToDevice(ligId, devices[i]);
		server.attachGridUnionToDevice(gridId, devices[i]);
		server.attachParamTableToDevice(devices[i]);
		server.addGPUWorker(devices[i]);
	}

	/* Finalize server initialization */
	if (devices.size() > 0) {
		server.updateDeviceIDLookup();
	}

	/*
	 ** The server and the data is now initialized. We are ready to use it.
	 ** Next we need to initialize the client.
	 */

	/* Allocate result buffer and declare dof buffer */
	vector<as::EnGrad> enGradOld(numDofs);
	as::EnGrad* enGradBuffer = enGradOld.data();

	as::DOF* dofBuffer = DOF_molecules[1].data();


	unsigned idx[2] = {0, 1};
	int reqIds[2] = {-1, -1};
	int reqIdsFirstIter[2] = {-1, -1}; // only used in first iteration
	unsigned DOFSize[2] = { (unsigned)ceil((double)numDofs/2.0), (unsigned)floor((double)numDofs/2.0) };

	assert(DOFSize[0]+DOFSize[1] == numDofs);
	assert(DOFSize[0] != 0);
	assert(DOFSize[1] != 0);
//	cout << "size0 " << DOFSize[0] << " size1 " <<  DOFSize[1] << endl;

	/* Devide initial DOF Buffer in two parts */
	vector<as::DOF> dof(numDofs);
	as::DOF* newDOF = dof.data();
	as::DOF* newDOFs[2] = {newDOF , newDOF + DOFSize[0]};
	as::DOF* oldDOFs[2] = {dofBuffer, dofBuffer + DOFSize[0]};

	vector<as::EnGrad> enGradNew(numDofs);
	as::EnGrad* newEnGrad = enGradNew.data();
	as::EnGrad* newEnGrads[2] = {newEnGrad, newEnGrad + DOFSize[0]};
	as::EnGrad* oldEnGrads[2] = {enGradBuffer, enGradBuffer + DOFSize[0]};



	/******** Main Loop ********/

	/* First iteration is unique since we cannot process anything meanwhile.
	 * The loop is of size 2 since we have two half buffers. */
	nvtxRangePushA("Processing");
	for (int i = 0; i < 2; ++i) {
		/* Submit Request for first half of DOF Buffer to the GPU Server.
		 * This is a non-blocking call */
		reqIdsFirstIter[idx[0]] = asClient::server_submit(server, oldDOFs[idx[0]], DOFSize[idx[0]], gridId, recId, ligId, serverMode);

		/* while energy is evaluated by the server, go ahead with calculating new configurations
		 * for first half of DOF Buffer */
		for(unsigned j = 0; j < DOFSize[idx[0]]; ++j) {
			const as::DOF& oldDOF = oldDOFs[idx[0]][j];
			as::DOF& newDOF = newDOFs[idx[0]][j];
			randomStep(oldDOF, newDOF);
		}

		/* Submit Request for first half of new DOF Buffer to the GPU Server.
		 * This is a non-blocking call */
		reqIds[idx[0]] = asClient::server_submit(server, newDOFs[idx[0]], DOFSize[idx[0]], gridId, recId, ligId, serverMode);

		/* Swap Buffers.
		 * The next use of idx[0] has the value of idx[1] and vice versa */
		std::swap(idx[0], idx[1]);
	}
	nvtxRangePop();

	/* Enter main loop */
	for(unsigned i = 0; i < (numIter-1)*2; ++i) {

		/* if we process a half buffer the first time, we need to wait for energies
		 * of the very first submission */
		nvtxRangePushA("Waiting");
		if (i == 0 || i == 1) {
			asClient::server_pull(server, reqIdsFirstIter[idx[0]], oldEnGrads[idx[0]]);
		}

		/* pull (wait) for request that was submitted two iterations ago
		 * the buffers */
		asClient::server_pull(server, reqIds[idx[0]], newEnGrads[idx[0]]);
		nvtxRangePop();

		/* Accept new positions according to Metropolis criterion */
		nvtxRangePushA("Processing");
		for(unsigned j = 0; j < DOFSize[idx[0]]; ++j) {
			as::DOF& oldDOF = oldDOFs[idx[0]][j];
			as::EnGrad& oldEnGrad = oldEnGrads[idx[0]][j];

			as::DOF& newDOF = newDOFs[idx[0]][j];
			as::EnGrad& newEnGrad = newEnGrads[idx[0]][j];

			/* the accepted values are stored in old variables !!! */
			MC_accept(oldDOF, oldEnGrad, newDOF, newEnGrad);

			if (idx[0] == 1 && j == 0 && false) {
				cout << "Iter "<< i/2 << endl;
				cout << oldDOF << endl;
				cout << oldEnGrad.E_El + oldEnGrad.E_VdW << endl;
			}

			/* calulate new trial configuration */
			randomStep(oldDOF, newDOF);
		}


		/* Submit Request to the GPU Server.
		 * This is a non-blocking call */
		reqIds[idx[0]] = asClient::server_submit(server, newDOFs[idx[0]], DOFSize[idx[0]], gridId, recId, ligId, serverMode);


		/* Swap Buffers */
		std::swap(idx[0], idx[1]);
		nvtxRangePop();

	}

	/* Finish last iteration */

	/* pull (wait) for request that was submitted two iterations ago
	 * the buffers */
	for (int i = 0; i < 2; ++i) {
		nvtxRangePushA("Waiting");
		if (numIter == 1) {
			asClient::server_pull(server, reqIdsFirstIter[idx[0]], oldEnGrads[idx[0]]);
		}

		asClient::server_pull(server, reqIds[idx[0]], newEnGrads[idx[0]]);
		nvtxRangePop();

		/* Accept new positions according to Metropolis criterion */
		nvtxRangePushA("Processing");
		for(unsigned j = 0; j < DOFSize[idx[0]]; ++j) {
			as::DOF& oldDOF = oldDOFs[idx[0]][j];
			as::EnGrad& oldEnGrad = oldEnGrads[idx[0]][j];

			as::DOF& newDOF = newDOFs[idx[0]][j];
			as::EnGrad& newEnGrad = newEnGrads[idx[0]][j];

			/* the accepted values are stored in old variables !!! */
			MC_accept(oldDOF, oldEnGrad, newDOF, newEnGrad);

		}
		nvtxRangePop();
		std::swap(idx[0], idx[1]);
	}


	/******** End Main Loop ********/


	/* print results to stderr */
	printResultsOutput(numDofs, dofBuffer, enGradBuffer, pivots);


	/* remove all data from host and devices */
	server.removeClient(clientId);


	log->info() << "exit" << endl;

	return 0;
}


