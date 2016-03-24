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
#include <fstream>
#include <nvToolsExt.h>

#include "as/ServerManagement.h"
#include "asUtils/Logger.h"

using namespace Log;

/*
 ** @brief: initializing the global Logger with file stream.
 */

void init_global_Logger(bool use_file = true) {
#ifndef NDEBUG
	logLevel level = Debug;
#else
	logLevel level = Info;
#endif

	if (use_file) {
		std::string filename = "AttractServer.log";
		global_log = new Logger(level, filename.substr(0,filename.size()-4));
	} else {
		global_log = new Logger(level, &(std::cerr));
	}

}

/* Constructor */
as::ServerManagement::ServerManagement() :
	W_mngt(*this, ITEMSIZE, DEVBUFSIZE, HOSTBUFSIZE), B_mngt(NUMITEMS, ITEMSIZE),
	Disp(*this, ITEMSIZE)
{
	init_global_Logger();
	global_log->info() << "ServerManagement.cpp: "<< "Starting Dispatcher" << std::endl;
	Disp.start();
}

as::ServerManagement::ServerManagement(unsigned itemSize) :
	W_mngt(*this, itemSize, DEVBUFSIZE, HOSTBUFSIZE), B_mngt(NUMITEMS, itemSize),
	Disp(*this, itemSize)
{
	init_global_Logger();
	global_log->info() << "ServerManagement.cpp: "<< "Starting Dispatcher" << std::endl;
	Disp.start();
}

as::ServerManagement::ServerManagement(unsigned numItems, unsigned itemSize) :
	W_mngt(*this, itemSize, DEVBUFSIZE, HOSTBUFSIZE), B_mngt(numItems, itemSize),
	Disp(*this, itemSize)
{
	init_global_Logger();
	global_log->info() << "ServerManagement.cpp: "<< "Starting Dispatcher" << std::endl;
	Disp.start();
}

as::ServerManagement::ServerManagement(unsigned numItems, unsigned itemSize, unsigned deviceBufferSize) :
	W_mngt(*this, itemSize, deviceBufferSize, HOSTBUFSIZE), B_mngt(numItems, itemSize),
	Disp(*this, itemSize)
{
	init_global_Logger();
	global_log->info() << "ServerManagement.cpp: "<< "Starting Dispatcher" << std::endl;
	Disp.start();
}

/* Destructor */
as::ServerManagement::~ServerManagement()
{
	global_log->info() << "ServerManagement.cpp: " << "Shutdown AttractServer: " << std::endl;
//	*global_log << std::setw(5) << " ";
	Disp.signalTerminate();
	Disp.join();
	global_log->info() << std::setw(5) << "Shutdown AttractServer: Ok" << std::endl;
}



/****************************
 * public member functions
 ****************************/
int as::ServerManagement::addProtein(int clientId, std::string filename)
{
	return D_mngt.addProtein(clientId, filename);
}

int as::ServerManagement::addProtein(std::string protein)
{
	return D_mngt.addProtein(protein);
}
int as::ServerManagement::addGridUnion(int clientId, std::string filename)
{
	return D_mngt.addGridUnion(clientId, filename);
}
int as::ServerManagement::addGridUnion(std::string gridUnion)
{
	return D_mngt.addGridUnion(gridUnion);
}

void as::ServerManagement::addParamTable(std::string name)
{
	D_mngt.addParamTable(name);
}

void as::ServerManagement::addSimParam(const SimParam& simPar) {
	D_mngt.addSimParam(simPar);
}

void as::ServerManagement::removeProtein(int globId)
{
	D_mngt.removeProtein(globId);
}

void as::ServerManagement::removeGridUnion(int globId)
{
	D_mngt.removeGridUnion(globId);
}

void as::ServerManagement::removeClient (int clientId)
{
	D_mngt.removeClient (clientId);
}

void as::ServerManagement::attachGridUnionToDevice(int globalId, unsigned deviceId)
{
	D_mngt.attachGridUnionToDevice(globalId, deviceId);
}

void as::ServerManagement::detachGridUnionFromDevice(int globalId, unsigned deviceId)
{
	D_mngt.detachGridUnionFromDevice(globalId, deviceId);
}

void as::ServerManagement::attachProteinToDevice(int globalId, unsigned deviceId)
{
	D_mngt.attachProteinToDevice(globalId, deviceId);
}

void as::ServerManagement::detachProteinFromDevice(int globalId, unsigned deviceId)
{
	D_mngt.detachProteinFromDevice(globalId, deviceId);
}

void as::ServerManagement::attachParamTableToDevice(unsigned deviceId)
{
	D_mngt.attachParamObjectsToDevice(deviceId);
}

void as::ServerManagement::detachParamTableFromDevice(unsigned deviceId)
{
	D_mngt.detachParamObjectsFromDevice(deviceId);
}

void as::ServerManagement::releaseDevice(unsigned deviceId)
{
	D_mngt.releaseDevice(deviceId);
}

void as::ServerManagement::updateDeviceIDLookup () {
	D_mngt.updateDeviceIDLookup();
}


std::string as::ServerManagement::dataInfo(dataInfoType_t infoType)
{
	return D_mngt.info(infoType);
}


void as::ServerManagement::addGPUWorker(unsigned deviceId)
{
	W_mngt.addGPUWorker(deviceId);
}

void as::ServerManagement::removeGPUWorker(unsigned deviceId)
{
	W_mngt.removeGPUWorker(deviceId);
}

void as::ServerManagement::addCPUWorker() {
	W_mngt.addCPUWorker();
}
void as::ServerManagement::removeCPUWorker() {
	W_mngt.removeCPUWorker();
}

int as::ServerManagement::submitRequest(DOF* dofs, const unsigned& numDOFs,
			const int& gridId,
			const int& recId, const int& ligId,
			const Request::useMode_t& mode)
{
	return Disp.createRequest(dofs, numDOFs,
			gridId, recId, ligId, mode);
}

as::Dispatcher::pullErr_t as::ServerManagement::pullRequest(int RequestId, EnGrad* buffer)
{
	return Disp.pullRequest(RequestId, buffer);
}



/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


