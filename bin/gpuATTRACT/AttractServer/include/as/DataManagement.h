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

#ifndef DATAMANAGEMENT_H_
#define DATAMANAGEMENT_H_

#include <string>
#include <vector>
#include <deque>
#include <set>
#include <algorithm>
#include <fstream>
#include <mutex>

#include "as/ClientData.h"
#include "as/Protein.h"
#include "as/GridUnion.h"
#include "as/ParamTable.h"
#include "as/SimParam.h"
#include "as/Shared.h"
#include "as/OccupancyLayout.h"
#include "asUtils/macros.h"


namespace as {

enum dataInfoType_t {
	client 	= 1,
	host 	= 2,
	device 	= 4,
	worker 	= 8
};

class mngt_io;

// TODO: make data manipulating functions threadsafe

class DataManagement {
public:
	/* Constructor */
	DataManagement(); // maybe reserve MAXCLIENT for vectors to avoid realloction
					 	// and therefore invalidation iterators. Realloction might lead to
						// invalid accesses from other threads.

	/* Destructor */
	~DataManagement();

	/* frieds */
	friend class mngt_io;

	/***************
	* G E T T E R
	***************/

	int globalProteinID(const int& clientId, const int& clientLocalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _clientData[clientId].globalProteinId(clientLocalId);
	}

	int globalGridID(const int& clientId, const int& clientLocalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _clientData[clientId].globalGridId(clientLocalId);
	}

	int deviceLocalProteinID(const int& deviceId, const int& globalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _prot_deviceOCC[deviceId].getLocation(globalId);
	}

	int deviceLocalGridID(const int& deviceId, const int& globalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _grid_deviceOCC[deviceId].getLocation(globalId);
	}

	const std::set<unsigned>& commonDeviceIDs(const int& globalGridId, const int& globalProteinId0, const int& globalProteinId1) {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		int minProtId = MIN(globalProteinId0, globalProteinId1);
		int maxProtId = MAX(globalProteinId0, globalProteinId1);
		const unsigned sizeSlice = (_prot_deviceIDs.size() * (_prot_deviceIDs.size() + 1)) / 2;
		const unsigned shift = (maxProtId*(maxProtId + 1)) / 2;
		std::lock_guard<std::mutex> guard(_m_devIDLookUp);
		return _commonDeviceID_LookUp[globalGridId*sizeSlice + shift + minProtId];
	}

	Protein* getProtein(const int& globalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _sharedProteins[globalId].get();
	}

	GridUnion* getGridUnion(const int& globalId) const {
		/* ToDo make runtime check in Debug mode: #ifndef NDEBUG ... if id is valid */
		return _sharedGridUnions[globalId].get();
	}

	AttrParamTable* getParamTable() const {
		return _sharedParamTable.get();
	}

	SimParam* getSimParam() const {
		return _sharedSimParam.get();
	}

	bool gridUnionValid(const int& globalId) {
		return ( ( globalId < static_cast<int>(_sharedGridUnions.size()) ) && _sharedGridUnions[globalId].isValid() );
	}

	bool proteinValid(const int& globalId) {
		return ( ( globalId < static_cast<int>(_sharedProteins.size()) ) && _sharedProteins[globalId].isValid() );
	}

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	/* return global id */
	int addProtein(int clientId, std::string filename);
	int addProtein (std::string protein);
	std::vector<int> addProteinEnsemble(int clientId, std::string filename);
	std::vector<int> addProteinEnsemble(std::string filename);

	int addGridUnion(int clientId, std::string filename);
	int addGridUnion(std::string gridUnion);

	void addParamTable(std::string name);
	void addSimParam(const SimParam& simPar);

	void removeProtein(int globalId);
	void removeGridUnion(int globalId);

	void removeClient (int clientId);

	void attachGridUnionToDevice(int globalId, unsigned deviceId);
	void detachGridUnionFromDevice(int globalId, unsigned deviceId);

	void attachProteinToDevice(int globalId, unsigned deviceId);
	void detachProteinFromDevice(int globalId, unsigned deviceId);

	/* Parameter Objects: parameter table (AttrParamTable)
	 * and simulation parameter (simParam) */
	void attachParamObjectsToDevice(unsigned deviceId);
	void detachParamObjectsFromDevice(unsigned deviceId);

	void releaseDevice(unsigned deviceId);

	std::string info(dataInfoType_t infoType);

	/*
	 ** @brief: This method has to be called after all proteins and grids are defined
	 */
	void updateDeviceIDLookup();

	/****************************
	 * public member variables
	 ****************************/

protected:
	/****************************
	 * protected member functions
	 ****************************/

	/****************************
	 * protected member variables
	 ****************************/

private:
	/****************************
	 * private member functions
	 ****************************/

	std::set<unsigned> clacCommonDeviceIDs(
			const int& globalGridId, const int& globalProtId0, const int& globalProtId1) const;

	std::set<unsigned> usedDevices();

	/****************************
	 * private member variables
	 ****************************/
	/** Contains global IDs of proteins and grids for each client */
	std::deque< ClientData > _clientData;

	/** Contain the actual grid/protein objects. They may be shared by multiple clients.
	 * The position in the container is the global ID of a grid/protein */
	std::deque< Shared<Protein> > _sharedProteins;
	std::deque< Shared<GridUnion> > _sharedGridUnions;
	Shared<AttrParamTable> _sharedParamTable;
	Shared<SimParam> _sharedSimParam;

	/** Contains a set of client IDs for each grid/protein.
	 * It tells which grid/protein is contained by which clients */
	std::deque< std::set<unsigned> > _grid_ClientIDs;
	std::deque< std::set<unsigned> > _prot_ClientIDs;

	/** Contains a set of devices IDs for each grid/protein/table/simParam.
	 * It tells which grid/protein/table/simParam is on which device */
	std::deque< std::set<unsigned> > _grid_deviceIDs;
	std::deque< std::set<unsigned> > _prot_deviceIDs;
	std::set<unsigned> _table_deviceIDs;
	std::set<unsigned> _simParam_deviceIDs;

	/** Lookup table for each Grid-Protein pair. The look-up returns a set of device IDs
	 * indicating the devices that can be used to process this pair.*/
	std::vector< std::set<unsigned> > _commonDeviceID_LookUp;

	/** Contains the occupancy layout for each device. */
	std::deque<OccupancyLayout> _grid_deviceOCC;
	std::deque<OccupancyLayout> _prot_deviceOCC;

	/** Contains the grid/protein device descriptions for each device */
	std::deque< std::vector<hostGridUnionResource> > _grid_deviceResc;
	std::deque< std::vector<hostProteinResource> > _prot_deviceResc;
	std::deque<hostParamTableResource> _table_deviceResc;

	/** mutexes */
	std::mutex _m_devIDLookUp;
	std::recursive_mutex _m_addData;

};

} // namespace



#endif // DATAMANAGEMENT_H_
