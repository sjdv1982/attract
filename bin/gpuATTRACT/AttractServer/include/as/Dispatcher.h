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

#ifndef DISPATCHER_H_
#define DISPATCHER_H_

#include <climits>
#include <atomic>
#include <list>
#include <unordered_map>
#include <queue>
#include <set>
#include <mutex>
#include <condition_variable>

#include "asUtils/Thread.h"
#include "as/asTypes.h"
#include "as/interface.h"
#include "as/Request.h"


namespace as {

class ServerManagement;

class Dispatcher : public asUtils::Thread {
public:
	/* Constructor */
	Dispatcher(ServerManagement& S_mngt, unsigned itemSize);

	/* Destructor */
	virtual ~Dispatcher();

	/* Types */
	enum pullErr_t {
		ready = 0,
		notReady = 1,
		invalid = 2,
	};

	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	void run();

	/*
	 ** @brief: creates Request from incomming DOF buffer. This function is
	 ** called by the client and is non-blocking despite of buffer
	 ** copying.
	 **
	 ** return: id of the request created. Can be used to query completeness.
	 ** -1 is returned if request cannot be accepted!
	 */
	int createRequest(DOF* dofs, const unsigned& numDOFs,
			const int& gridId, const int& recId, const int& ligId,
			const Request::useMode_t& mode);

	pullErr_t pullRequest(const int& reqId, EnGrad* const buffer);

	/*
	 ** @brief: Signal the Dispatcher to terminate. It does not terminate as long as
	 ** there are requests in the queue;
	 */
	inline void signalTerminate() {
		_terminate = true;
		_condVar.notify_one();
	}

	/*
	 ** @brief: this is done by the dispatcher itself
	 */
	inline Request* removeRequest();

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
	inline bool reqIdValid (const int& id) {
		/* request with key id is constained in the map */
		return _reqMap.find(id) != _reqMap.end();
	}

	inline WorkerItem* getItemFromBufMngt();

	inline void returnItemToBufMngt(WorkerItem* const item);

	inline bool gridValid(const int& globalId);

	inline bool protValid(const int& globalId);

	inline const std::set<unsigned>& getCommonDeviceIDs_cref (const int& globalGridId, const int& globalProteinId0, const int& globalProteinId1) const;

	inline unsigned deviceFillLevel(const unsigned deviceId) const;
	inline unsigned hostFillLevel(const unsigned id) const;

	inline void pushItemToGPUWorker(WorkerItem* const item, const unsigned& deviceId);

	/* Rule for Distributing the items across workers */
	inline void GPUDist(Request* const req);

	/* Rule for Distributing the items across workers */
	inline void CPUDist(Request* const req);

	inline void pushItemToCPUWorker(WorkerItem* const item, const unsigned& id);

	inline unsigned numCPUWorker() const;

	/****************************
	 * private member variables
	 ****************************/
	unsigned _itemSize;

	std::unordered_map<int, Request*> _reqMap;
	std::queue<Request*, std::list<Request*> > _reqQueue;

	std::atomic<bool> _terminate;
	std::atomic<bool> _sleeping;

	std::mutex _mutexMap;
	std::mutex _mutexQueue;
	std::condition_variable _condVar;

	std::atomic<int> _id;

	ServerManagement& _S_mngt;
};

} // namespace



#endif /* DISPATCHER_H_ */
