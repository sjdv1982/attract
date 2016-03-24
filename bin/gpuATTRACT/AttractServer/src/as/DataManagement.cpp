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
#include <iomanip>
#include <iostream>
#include <sstream>

#include <cstdlib>
#include <cassert>
#include <cmath>

#include "as/DataManagement.h"
#include "as/DeviceDataInterface.h"
#include "as/GridUnion.h"
#include "as/NLGrid.h"
#include "as/IntrplGrid.h"
#include "as/Protein.h"
#include "as/DeviceDataFactory.h"
#include "as/ParamTable.h"
#include "as/interface.h"
#include "as/mngt_io.h"
#include "asDB/readFile.h"
#include "config.h"
#include "asUtils/Logger.h"

using namespace Log;

/* Constructor */
as::DataManagement::DataManagement() :
	_sharedSimParam(new SimParam) {}

/* Destructor */
as::DataManagement::~DataManagement()
{
	global_log->info() << "Shutdown DataManagement:" << std::endl;
	/* Detach all device resources */
	/* The host resources get automatically freed by shared pointers */
	std::set<unsigned> devicesInUse = usedDevices();
	global_log->info() << std::setw(5) << " "<< "Freeing Resources on Device: ";
	for(auto i : devicesInUse) {
		*global_log << i << " ";
	}

	for(auto i : devicesInUse) {
		releaseDevice(i);
	}

	*global_log << std::endl;
	global_log->info() << std::setw(5) << " " << "Ok" << std::endl;

	/* Bad!!! but necessary. it is constructed in server management */
	delete global_log;
}

/****************************
 * public member functions
 ****************************/

/*
 ** @brief:
 */
int as::DataManagement::addProtein (int clientId, std::string name)
{
	/* search in global grid/protein arrays if already used and assign a global ID */
	/* clientId == -1 is a special case: It is used to add a Grid without holding it by a client */
	ASSERT(name.empty() == false);
	ASSERT(clientId >= -1);
	int globId = -1;
	{
		std::lock_guard<std::recursive_mutex> guard(_m_addData);
		if (clientId >= static_cast<int>(_clientData.size()) && clientId != -1) {
			assert(clientId == static_cast<int>(_clientData.size()));
			_clientData.push_back(ClientData());
		}

		/* check if Protein is already in data structure */
		for (int i = 0; i < static_cast<int>(_sharedProteins.size()); ++i) {
			if (_sharedProteins[i].isValid()) {
				if (_sharedProteins[i].get()->tag() == name) {
					globId = i;
					if (clientId != -1) {
						_sharedProteins[i].incrCount();
					}
					break;
				}
			}
		}
	}
	/* if not, read it and assign it to the first available location */
	if (globId == -1) {
		Protein* obj = asDB::createProteinFromPDB(name);
		/* auto-pivotize:
		 * TODO move that elsewhere. This is behavior that is not expected */
		obj->auto_pivotize();

		Shared<Protein> sharedObj(obj);

		std::lock_guard<std::recursive_mutex> guard(_m_addData);

		/* get the next free ID and assign protein */
		if (_sharedProteins.size() == 0) {
			_sharedProteins.push_back(sharedObj);
			globId = 0;
		} else {
			int i;
			for (i = 0; i < static_cast<int>(_sharedProteins.size()); ++i) {
				if (!_sharedProteins[i].isValid()) {
					_sharedProteins[i] = sharedObj;
					globId = i;
					break;
				}
			}
			if (globId == -1) {
				_sharedProteins.push_back(sharedObj);
				globId = _sharedProteins.size()-1;
			}
		}
	}

	assert(globId != -1);
	if (clientId != -1) {
		std::lock_guard<std::recursive_mutex> guard(_m_addData);
		_clientData[clientId].addProtein(globId);
		if (globId >= static_cast<int>(_prot_ClientIDs.size())) {
			assert(globId == static_cast<int>(_prot_ClientIDs.size()));
			_prot_ClientIDs.push_back(std::set<unsigned>());
			_prot_ClientIDs[globId].insert(clientId);
		} else {
			_prot_ClientIDs[globId].insert(clientId);
		}
	}

	return globId;
}

int as::DataManagement::addGridUnion (int clientId, std::string name) {
	/* search in global grid/protein arrays if already used and assign a global ID */
	/* clientId == -1 is a special case: It is used to add a Grid without holding it by a client */

	ASSERT(name.empty() == false);
	ASSERT(clientId >= -1);
	int globId = -1;
	{
		std::lock_guard<std::recursive_mutex> guard(_m_addData);
		if (clientId >= static_cast<int>(_clientData.size()) && clientId != -1) {
			assert(clientId == static_cast<int>(_clientData.size()));
			_clientData.push_back(ClientData());
		}

		/* check if grid exists already in data structure */
		for (unsigned i = 0; i < _sharedGridUnions.size(); ++i) {
			if (_sharedGridUnions[i].isValid()) {
				if (_sharedGridUnions[i].get()->tag() == name) {
					globId = i;
					if (clientId != -1) {
						_sharedGridUnions[i].incrCount();
					}
					break;
				}
			}
		}
	}
	/* if not, read it and assign it to the first available location */
	if (globId == -1) {
		GridUnion* obj = asDB::createGridFromGridFile(name);
		Shared<GridUnion> sharedObj(obj);

		std::lock_guard<std::recursive_mutex> guard(_m_addData);

		/* get the next free ID and assign gridUnion */
		if (_sharedGridUnions.size() == 0) {
			_sharedGridUnions.push_back(sharedObj);
			globId = 0;
		} else {
			int i;
			for (i = 0; i < static_cast<int>(_sharedGridUnions.size()); ++i) {
				if (!_sharedGridUnions[i].isValid()) {
					_sharedGridUnions[i] = sharedObj;
					globId = i;
					break;
				}
			}
			if (globId == -1) {
				_sharedGridUnions.push_back(sharedObj);
				globId = _sharedGridUnions.size()-1;
			}
		}
	}

	assert(globId != -1);
	if (clientId != -1) {
		std::lock_guard<std::recursive_mutex> guard(_m_addData);
		_clientData[clientId].addGrid(globId);

		if (globId >= static_cast<int>(_grid_ClientIDs.size())) {
			assert(globId == static_cast<int>(_grid_ClientIDs.size()));
			_grid_ClientIDs.push_back(std::set<unsigned>());
			_grid_ClientIDs[globId].insert(clientId);
		} else {
			_grid_ClientIDs[globId].insert(clientId);
		}
	}

	return globId;
}

int as::DataManagement::addProtein (std::string name) {
	return addProtein(-1, name);
}

std::vector<int> as::DataManagement::addProteinEnsemble(int clientId, std::string filename) {
	std::vector<std::string> fileNames = asDB::readFileNamesFromEnsembleList(filename);
	std::vector<int> proteinIds;
	for (auto i = 0; i < fileNames.size(); ++i) {
		int id = addProtein(clientId, fileNames[i]);
		proteinIds.push_back(id);
	}
	/* Pivotize explicitly since all need to have the same pivot*/
	Protein* prot0 = getProtein(proteinIds[0]);
	for(auto id: proteinIds) {
		getProtein(id)->pivotize(prot0->pivot());
	}
	return proteinIds;
}
std::vector<int> as::DataManagement::addProteinEnsemble(std::string filename) {
	return addProteinEnsemble(-1, filename);
}

int as::DataManagement::addGridUnion (std::string filename) {
	return addGridUnion(-1, filename);
}

void as::DataManagement::addParamTable(std::string name)
{
	ASSERT(name.empty() == false);
	AttrParamTable* obj = asDB::createParamTableFromFile(name);
	Shared<AttrParamTable> sharedObj(obj);

	std::lock_guard<std::recursive_mutex> guard(_m_addData);
	_sharedParamTable = sharedObj;

}

/*
 ** @brief: Adds a custom SimPar object.
 ** The default one is discarded by assignment operation
 ** of share object.
 */
void as::DataManagement::addSimParam(const SimParam& simPar) {
	SimParam* simPar_ptr = new SimParam(simPar);
	Shared<SimParam> sharedObj(simPar_ptr);
	std::lock_guard<std::recursive_mutex> guard(_m_addData);
	_sharedSimParam = sharedObj;
}

void as::DataManagement::removeProtein(int globId)
{
	/* check if global Id is valid */
	if ( (globId >= static_cast<int>(_sharedProteins.size()) )  || !_sharedProteins[globId].isValid()) {

		std::stringstream stream;
		stream << "Removing protein with global ID " << globId
				<< " failed.\n";
		stream << "Protein ID is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* destroys the object */
	_sharedProteins[globId].decrCount();

	/* remove grid from devices if necessary */
	/* assert that the grid is now invalid */
	assert(!_sharedProteins[globId].isValid());

	if ( globId < static_cast<int>(_prot_deviceIDs.size()) ) {
		const std::set<unsigned> devices = _prot_deviceIDs[globId];
		for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
			detachProteinFromDevice(globId, *it);
		}
	}
}

/*
 ** @brief:
 ** Convention: The grid gets automatically removed also from all devices.
 */
void as::DataManagement::removeGridUnion(int globId)
{
	/* check if global Id is valid */
	if ( (globId >= static_cast<int>(_sharedGridUnions.size()))  || !_sharedGridUnions[globId].isValid()) {
		std::stringstream stream;
		stream << "Removing grid with global ID " << globId
				<< " failed.\n";
		stream << "Grid ID is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* destroys the object */
	_sharedGridUnions[globId].decrCount();

	/* remove grid from devices if necessary */
	/* assert that the grid is now invalid */
	assert(!_sharedGridUnions[globId].isValid());

	if ( globId < static_cast<int>(_grid_deviceIDs.size())) {
		const std::set<unsigned> devices = _grid_deviceIDs[globId];
		for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
			detachGridUnionFromDevice(globId, *it);
		}
	}

}


void as::DataManagement::removeClient (int clientId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if clientId is valid */
	if(clientId >= static_cast<int>(_clientData.size()) || !_clientData[clientId].isValid()) {
		std::stringstream stream;
		stream << "Removing Client is not possible. Client ID is not valid" << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	ClientData& client = _clientData[clientId];
	std::vector<int> proteins = client.proteins();
	for (unsigned i = 0; i<proteins.size(); ++i) {
		int globId = proteins[i];
		_sharedProteins[globId].decrCount();
		_prot_ClientIDs[globId].erase(clientId);

		/* remove protein from devices if necessary */
		if (!_sharedProteins[globId].isValid() && (globId < static_cast<int>(_prot_deviceIDs.size()) ) ) {
			// do not use a reference because the elements of the set are removed inside detachGridUnionFromDevice
			std::set<unsigned> devices = _prot_deviceIDs[globId];
			for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
				detachProteinFromDevice(globId, *it);
			}
		}
	}
	std::vector<int> grids = client.grids();
	for (unsigned i = 0; i<grids.size(); ++i) {
		int globId = grids[i];
		_sharedGridUnions[globId].decrCount();
		_grid_ClientIDs[globId].erase(clientId);

		/* remove grid from devices if necessary */
		if (!_sharedGridUnions[globId].isValid() && (globId < static_cast<int>(_grid_deviceIDs.size())) ) {
			// do not use a reference because the elements of the set are removed inside detachGridUnionFromDevice
			const std::set<unsigned> devices = _grid_deviceIDs[globId];
			for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
				detachGridUnionFromDevice(globId, *it);
			}
		}
	}
	_clientData[clientId].reset();
}

void as::DataManagement::attachGridUnionToDevice(int globalId, unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if grid is valid */
	if ( (globalId >= static_cast<int>(_sharedGridUnions.size()))  || !_sharedGridUnions[globalId].isValid()) {
		std::stringstream stream;
		stream << "Attaching Grid with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Grid ID is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* check if device Id is valid */
	// ToDo

	if (deviceId >= _grid_deviceOCC.size()) {
		_grid_deviceOCC.resize(deviceId + 1, OccupancyLayout(DEVICE_MAXGRIDS));
		_grid_deviceResc.resize(deviceId + 1, std::vector<hostGridUnionResource>(DEVICE_MAXGRIDS));
	} else if (_grid_deviceOCC[deviceId].contains(globalId)) {
		/* check if grid is already attached to the device */
		std::stringstream stream;
		stream << "Attaching Grid with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Grid is already attached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* check if device limit is reached */
	if(_grid_deviceOCC[deviceId].isfull()) {
		std::stringstream stream;
		stream << "Attaching Grid with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Device maximum has been reached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	if (globalId >= static_cast<int>(_grid_deviceIDs.size())) {
		_grid_deviceIDs.resize(globalId + 1);
	}

	/* adapt the grid's deviceId list */
	_grid_deviceIDs[globalId].insert(deviceId);

	/* adapt device occupancy layout */
	int devloc = _grid_deviceOCC[deviceId].getFirstEmptyLocation();
	assert(devloc >= 0);
	_grid_deviceOCC[deviceId].addObject(globalId, devloc);

	GridUnion* obj = _sharedGridUnions[globalId].get();
	cudaGridUnionDesc cudaDesc = DeviceDataFactory::initDeviceGridUnion(obj, deviceId);
	deviceGridUnionDesc &deviceDesc = cudaDesc.deviceDesc;
	hostGridUnionResource &hostResc = cudaDesc.hostResc;
	setDeviceGridUnion(deviceDesc, deviceId, devloc);

	/* adapt device descriptions */
	std::vector<hostGridUnionResource>& vec = _grid_deviceResc[deviceId];
	vec[devloc] = hostResc;
}

void as::DataManagement::detachGridUnionFromDevice(int globalId, unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if device Id is valid */
	// ToDo

	/* check if protein is located at the device */
	if ((globalId >= static_cast<int>(_grid_deviceIDs.size())) ||
			(_grid_deviceIDs[globalId].find(deviceId) == _grid_deviceIDs[globalId].end()) ) {
		/* it is not contained ... */
		std::stringstream stream;
		stream << "Cannot detach Grid with global ID " << globalId << " from device "
				<< deviceId << " because it has not been attached" << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	int devloc = _grid_deviceOCC[deviceId].getLocation(globalId);

	std::vector<hostGridUnionResource>& vec = _grid_deviceResc[deviceId];
	hostGridUnionResource hostResc = vec[devloc];

	unsetDeviceGridUnion(deviceId, devloc);
	DeviceDataFactory::disposeDeviceGridUnion(hostResc, deviceId);

	_grid_deviceIDs[globalId].erase(deviceId);
	_grid_deviceOCC[deviceId].removeObj(globalId);
	vec[devloc] = hostGridUnionResource();
}

void as::DataManagement::attachProteinToDevice(int globalId, unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if protein is valid */
	if ( (globalId >= static_cast<int>(_sharedProteins.size())) || !_sharedProteins[globalId].isValid() ) {
		std::stringstream stream;
		stream << "Attaching Protein with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Protein ID is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* check if device Id is valid */
	// ToDo

	if (deviceId >= _prot_deviceOCC.size()) {
		_prot_deviceOCC.resize(deviceId + 1, OccupancyLayout(DEVICE_MAXPROTEINS));
		_prot_deviceResc.resize(deviceId + 1, std::vector<hostProteinResource>(DEVICE_MAXPROTEINS));
	} else if (_prot_deviceOCC[deviceId].contains(globalId)) {
		/* check if grid is already attached to the device */
		std::stringstream stream;
		stream << "Attaching Protein with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Protein is already attached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* check if device limit is reached */
	if(_prot_deviceOCC[deviceId].isfull()) {
		std::stringstream stream;
		stream << "Attaching Protein with global ID " << globalId
				<< " to Device " << deviceId << " failed.\n";
		stream << "Device maximum has been reached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	if (globalId >= static_cast<int>(_prot_deviceIDs.size())) {
		_prot_deviceIDs.resize(globalId + 1);
	}

	/* adapt the protein's deviceId list */
	_prot_deviceIDs[globalId].insert(deviceId);

	/* adapt device occupancy layout */
	int devloc = _prot_deviceOCC[deviceId].getFirstEmptyLocation();
	assert(devloc >= 0);
	_prot_deviceOCC[deviceId].addObject(globalId, devloc);

	Protein* obj = _sharedProteins[globalId].get();
	cudaProteinDesc cudaDesc = DeviceDataFactory::initDeviceProtein(obj, deviceId);

	deviceProteinDesc &deviceDesc = cudaDesc.deviceDesc;
	hostProteinResource &hostResc = cudaDesc.hostResc;

	setDeviceProtein(deviceDesc, deviceId, devloc);

	/* adapt device descriptions */
	std::vector<hostProteinResource>& vec = _prot_deviceResc[deviceId];
	vec[devloc] = hostResc;
}

void as::DataManagement::detachProteinFromDevice(int globalId, unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if device Id is valid */
	// ToDo

	/* check if protein is located at the device */
	if ((globalId >= static_cast<int>(_prot_deviceIDs.size()) ) ||
			(_prot_deviceIDs[globalId].find(deviceId) == _prot_deviceIDs[globalId].end()) ) {
		/* it is not contained ... */
		std::stringstream stream;
		stream << "Cannot detach Protein with global ID " << globalId << " from device "
				<< deviceId << " because it has not been attached" << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	int devloc = _prot_deviceOCC[deviceId].getLocation(globalId);

	std::vector<hostProteinResource>& vec = _prot_deviceResc[deviceId];
	hostProteinResource hostResc = vec[devloc];

	unsetDeviceProtein(deviceId, devloc);
	DeviceDataFactory::disposeDeviceProtein(hostResc, deviceId);

	_prot_deviceIDs[globalId].erase(deviceId);
	_prot_deviceOCC[deviceId].removeObj(globalId);
	vec[devloc] = hostProteinResource();
}

void as::DataManagement::attachParamObjectsToDevice(unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if device Id is valid */
	// ToDo

	/* attach table */

	/* check if table is already attached to the device */
	if (_table_deviceIDs.find(deviceId) != _table_deviceIDs.end()) {
		std::stringstream stream;
		stream << "Attaching table to Device " << deviceId << " failed.\n";
		stream << "Table is already attached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	if (_sharedParamTable.isValid()) {
		_table_deviceIDs.insert(deviceId);
		AttrParamTable* obj = _sharedParamTable.get();
		cudaParamTableDesc cudaDesc = DeviceDataFactory::initDeviceParamTable(obj, deviceId);
		deviceParamTableDesc &deviceDesc = cudaDesc.deviceDesc;
		hostParamTableResource &hostResc = cudaDesc.hostResc;
		setDeviceParamTable(deviceDesc, deviceId);
		if (deviceId >= _table_deviceResc.size()) {
			_table_deviceResc.resize(deviceId + 1);
		}
		_table_deviceResc[deviceId] = hostResc;
	} else {
		std::stringstream stream;
		stream << "Attaching table with to Device " << deviceId << " failed.\n";
		stream << "Table is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	/* attach SimParam */

	/* check if SimParam is already attached to the device */
	if (_simParam_deviceIDs.find(deviceId) != _simParam_deviceIDs.end()) {
		std::stringstream stream;
		stream << "Attaching SimParam to Device " << deviceId << " failed.\n";
		stream << "SimParam is already attached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	if (_sharedSimParam.isValid()) {
		_simParam_deviceIDs.insert(deviceId);
		SimParam* obj = _sharedSimParam.get();
		/* SimParam is of POD type -> usable at both the device and host -> no host/device-Desc type neccessary */

		setDeviceSimParam(*obj, deviceId);
	} else {
		std::stringstream stream;
		stream << "Attaching SimParam to Device " << deviceId << " failed.\n";
		stream << "SimParam is not valid." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}
}
void as::DataManagement::detachParamObjectsFromDevice(unsigned deviceId)
{
	std::lock_guard<std::recursive_mutex> guard(_m_addData);

	/* check if device Id is valid */
	// ToDo

	/* check if table is located at the device */
	if (_table_deviceIDs.find(deviceId) == _table_deviceIDs.end()) {
		std::stringstream stream;
		stream << "Detaching table from to Device " << deviceId << " failed.\n";
		stream << "Table has not been attached." << std::endl;
		global_log->error() << stream.str() << std::endl;
		return;
	}

	hostParamTableResource hostResc = _table_deviceResc[deviceId];
	_table_deviceResc[deviceId] = hostParamTableResource();

	DeviceDataFactory::disposeDeviceParamTable(hostResc, deviceId);
	_table_deviceIDs.erase(deviceId);

}

void as::DataManagement::updateDeviceIDLookup() {
	const int sizeGrids = _grid_deviceIDs.size();
	const int sizeProts = _prot_deviceIDs.size();
	const int sizeSlice = (sizeProts * (sizeProts + 1)) / 2;
	const int size =  sizeGrids * sizeSlice;

	/* lookup layout: lower trianglular matrix */

	std::vector< std::set<unsigned> > lookUp(size);
	for (int i = 0; i < sizeGrids; ++i) {
		for (int j = 0; j < sizeProts; ++j) {
			const int shift = (j*(j+1)) / 2;
			for(int k = 0; k <= j; ++k) {
				lookUp[i*sizeSlice + shift + k] = clacCommonDeviceIDs(i,j,k);

			}
		}
	}
	std::lock_guard<std::mutex> guard(_m_devIDLookUp);
	_commonDeviceID_LookUp = lookUp;
}

void as::DataManagement::releaseDevice(unsigned deviceId) {
	/* Remove all Proteins, Grids and Tables */
	std::lock_guard<std::recursive_mutex> guard(_m_addData);
	std::set<int> globalIds;
	globalIds = _prot_deviceOCC[deviceId].objectSet();

//	std::cout << "DataManagement.cpp: " << "global_protId on Device: " << deviceId << " : ";
//	for (auto globId : globalIds) {
//		std::cout << globId  << " ";
//	}
//	std::cout << std::endl;
	for (int globId : globalIds) {
		detachProteinFromDevice(globId, deviceId);
	}

	globalIds = _grid_deviceOCC[deviceId].objectSet();
	for (int globId : globalIds) {
		detachGridUnionFromDevice(globId, deviceId);
	}

	detachParamObjectsFromDevice(deviceId);
}


std::string as::DataManagement::info(dataInfoType_t infoType)
{
	/* TODO: implement host, device, worker output */
	std::string str;

	switch (infoType) {
	case client:
		str = mngt_io::clientInfo(this);
		break;
	case host:
		break;
	case device:
		break;
	case worker:
		break;
	}
	return str;
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


std::set<unsigned> as::DataManagement::clacCommonDeviceIDs(const int& globalGridId,
		const int& globalProtId0, const int& globalProtId1) const {
	const std::set<unsigned>& grid = _grid_deviceIDs[globalGridId]; 	/** grid device IDs*/
	const std::set<unsigned>& prot0 = _prot_deviceIDs[globalProtId0]; 	/** protein device IDs*/
	const std::set<unsigned>& prot1 = _prot_deviceIDs[globalProtId1];
	std::set<unsigned> intersect0;
	std::set<unsigned> intersect1;
	std::set_intersection(grid.begin(), grid.end(), prot0.begin(), prot0.end(),
			std::insert_iterator< std::set<unsigned> >(intersect0, intersect0.begin())
			);
	std::set_intersection(intersect0.begin(), intersect0.end(), prot1.begin(), prot1.end(),
			std::insert_iterator< std::set<unsigned> >(intersect1, intersect1.begin())
			);
	return intersect1;
}


std::set<unsigned> as::DataManagement::usedDevices() {

	/* built unions of ..._deviceIds by inserting all sets
	 * duplicated are removed because std::set contains only
	 * unique elements */
	std::set<unsigned> result;
	for (auto set : _prot_deviceIDs) {
		result.insert(set.begin(), set.end());
	}
	for (auto set : _grid_deviceIDs) {
		result.insert(set.begin(), set.end());
	}
	result.insert(_table_deviceIDs.begin(), _table_deviceIDs.end());
	return result;
}
