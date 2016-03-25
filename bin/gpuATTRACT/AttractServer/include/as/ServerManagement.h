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

#ifndef SERVERMANAGEMENT_H_
#define SERVERMANAGEMENT_H_

#include "as/DataManagement.h"
#include "as/WorkerManagement.h"
#include "as/BufferManagement.h"
#include "as/Dispatcher.h"
#include "as/SimParam.h"
#include "as/interface.h"
#include "as/Request.h"
#include "config.h"

namespace as {


class ServerManagement {
public:
	/* Constructor */
	ServerManagement();
	ServerManagement(unsigned itemSize);
	ServerManagement(unsigned numItems, unsigned itemSize);
	ServerManagement(unsigned numItems, unsigned itemSize, unsigned deviceBufferSize);

	/* Destructor */
	~ServerManagement();

	/***************
	* G E T T E R
	***************/

	DataManagement* DataMngt() {
		return &D_mngt;
	}

	BufferManagement* BufferMngt() {
		return &B_mngt;
	}

	WorkerManagement* WorkerMngt() {
		return &W_mngt;
	}

	Protein* getProtein(const int& globalId) const {
		return D_mngt.getProtein(globalId);
	}

	GridUnion* getGridUnion(const int& globalId) const {
		return D_mngt.getGridUnion(globalId);
	}

	AttrParamTable* getParamTable() const {
		return D_mngt.getParamTable();
	}

	SimParam* getSimParam() const {
		return D_mngt.getSimParam();
	}

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/

	int addProtein (int clientId, std::string filename);
	int addProtein (std::string protein);
	int addGridUnion (int clientId, std::string filename);
	int addGridUnion (std::string gridUnion);

	void addParamTable(std::string name);
	void addSimParam(const SimParam& simPar);

	void removeProtein(int globId);
	void removeGridUnion(int globId);

	void removeClient (int clientId);

	void attachGridUnionToDevice(int globalId, unsigned deviceId);
	void detachGridUnionFromDevice(int globalId, unsigned deviceId);

	void attachProteinToDevice(int globalId, unsigned deviceId);
	void detachProteinFromDevice(int globalId, unsigned deviceId);

	void attachParamTableToDevice(unsigned deviceId);
	void detachParamTableFromDevice(unsigned deviceId);

	void updateDeviceIDLookup();
	void releaseDevice(unsigned deviceId);

	std::string dataInfo(dataInfoType_t infoType);

	void addGPUWorker(unsigned deviceId);
	void removeGPUWorker(unsigned deviceId);

	void addCPUWorker();
	void removeCPUWorker();

	int submitRequest(DOF* dofs, const unsigned& numDOFs,
			const int& gridId, const int& recId, const int& ligId,
			const Request::useMode_t& mode);

	Dispatcher::pullErr_t pullRequest(int RequestId, EnGrad* buffer);


	/* TODO: Implementation:
	 * Someday a GPU Worker shall dynamically be removed, the corresponding
	 * device id has to be removed from respective containers.
	 * (_grid_deviceIDs, _prot_deviceIDs)
	 * After that, the lookup table for common device ids has to be updated
	 * so that no workerItems are distributed to that worker anymore.
	 * The Worker can safely be terminated after having processed the remaining work.
	 * After the worker terminated, all Grids, Proteins and the table can be detached
	 * from the device.
	 * */
	void releaseGPU(unsigned deviceId);


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

	/****************************
	 * private member variables
	 ****************************/
	DataManagement D_mngt;
	WorkerManagement W_mngt;
	BufferManagement B_mngt;
	Dispatcher Disp;

};

} // namespace

#endif /* SERVERMANAGEMENT_H_ */
