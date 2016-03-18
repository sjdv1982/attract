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

#include <cstring>
#include <cassert>
#include <set>
#include <algorithm>
#include <iostream>

#include "as/Dispatcher.h"
#include "as/ServerManagement.h"
#include "asUtils/Logger.h"

#include "config.h"

using namespace Log;


/* Constructor */
as::Dispatcher::Dispatcher(ServerManagement& S_mngt, unsigned itemSize) :
		_itemSize(itemSize), _terminate(false), _sleeping(false), _id(0), _S_mngt(S_mngt) {}


/* Destructor */
as::Dispatcher::~Dispatcher() {
	global_log->info() << "Shutdown Dispatcher: " << std::endl;
	if(_sleeping == true || _terminate == false) {
		global_log->info() << std::setw(5) << " " << "Error" << std::endl;
	} else {
		global_log->info() << std::setw(5) << " " << "Ok" << std::endl;
	}
}


/****************************
 * public member functions
 ****************************/
void as::Dispatcher::run() {
	while (true) {
		/* fetch a request */
		Request* req;
		req = removeRequest(); /* goes to sleep if no request is available */
		if (req == NULL) {
			/* It was signaled that dispatcher has to terminate */
			break;
		}
		switch(req->useMode()) {
		case Request::GPU:
			/* Analyze the filling level of each the WorkQueue */
			GPUDist(req);
			break;
		case Request::CPU:
			/* Distribute items among CPU workers */
			CPUDist(req);
			break;
		case Request::GPUCPU:
			/* Distribute items among CPU and GPU workers */
			global_log->error() << "Hybrid GPU-CPU use-mode not supported" << std::endl;
			std::exit(EXIT_FAILURE);
			break;
		case Request::unspecified:
			/* intentionally left blank */
		default:
			global_log->error() << "Use-Mode unspecified" << std::endl;
			std::exit(EXIT_FAILURE);

		}

	}
	global_log->info() << "Dispatcher terminated" << std::endl;
}

inline void as::Dispatcher::GPUDist(Request* const req)
{
	const unsigned& globGridId = req->globGridId();
	const unsigned& globRecId = req->globRecId();
	const unsigned& globLigId = req->globLigId();
	const std::set<unsigned>& commonDevices = getCommonDeviceIDs_cref(globGridId, globRecId, globLigId);

	const unsigned numDevices = commonDevices.size();
	assert(numDevices > 0);

	/* reset the iterator to call req->nextItem. Accordingly, it starts to return the first item*/
	req->resetIterator();

	if (numDevices == 1) {
		/* Push items directly to WorkQueue of the respective GPU worker */
		for (unsigned i = 0; i < req->size(); ++i) {
			pushItemToGPUWorker(req->nextItem(), *commonDevices.begin());
		}
	} else {
//		std::cout << "Dispatcher.cpp :" << "numDevices " << numDevices <<  std::endl;

		/* Analyze the filling level of each the WorkQueue */
		DevFillLevel devFillLevels[numDevices];


		unsigned count = 0;
		for (unsigned i : commonDevices) {
			devFillLevels[count].id = i;
			devFillLevels[count].level = deviceFillLevel(i);
//			std::cout << "Dispatcher.cpp :" << "id deviceFillLevel "
//					<< devFillLevels[count].id << devFillLevels[count].level <<  std::endl;
			++count;

		}

		/* Sort by increasing fill level */
		std::sort(devFillLevels, devFillLevels + numDevices);

//		count = 0;
//		std::cout << "Dispatcher.cpp :" << "numItems " << req->size() <<  std::endl;
//		std::cout << "Dispatcher.cpp :" << "Push Item " << count++
//							<<  " to Device " << devFillLevels[0].id <<  std::endl;

		/* push first item to Worker with the lowest level first */
		pushItemToGPUWorker(req->nextItem(), devFillLevels[0].id);

		/* The following statements conserve that devFillLevels stays sorted */
		for (unsigned i = 0; i < req->size()-1; ++i) {
			for (unsigned j = 0; j < numDevices-1; ++j) {
				if(devFillLevels[j] < devFillLevels[j + 1]) {

//					std::cout << "Dispatcher.cpp :" << "0 Push Item " << count++
//							<<  " to Device " << devFillLevels[j].id <<  std::endl;

					pushItemToGPUWorker(req->nextItem(), devFillLevels[j].id);
					++devFillLevels[j].level;



					break;
				}
				if (j == numDevices-2) {

//					std::cout << "Dispatcher.cpp :" << "1 Push Item " << count++
//							<<  " to Device " << devFillLevels[numDevices-1].id <<  std::endl;

					pushItemToGPUWorker(req->nextItem(), devFillLevels[numDevices-1].id);
					++devFillLevels[numDevices-1].level;
				}
			}
		}
	}

}

inline void as::Dispatcher::CPUDist(Request* const req)
{
	const unsigned numWorkers = numCPUWorker();
	/* reset the iterator to call req->nextItem. Accordingly, it starts to return the first item*/
	req->resetIterator();
	if (numWorkers == 1) {
		/* Push items directly to WorkQueue of the respective GPU worker */
		for (unsigned i = 0; i < req->size(); ++i) {
			WorkerItem* item = req->nextItem();
			pushItemToCPUWorker(item, 0);
		}
	} else {

		/* Analyze the filling level of each the WorkQueue */
		DevFillLevel fillLevels[numWorkers];


		unsigned count = 0;
		for (unsigned i = 0; i < numWorkers; ++i) {
			fillLevels[count].id = i;
			fillLevels[count].level = hostFillLevel(i);
//			std::cout << "Dispatcher.cpp :" << "id deviceFillLevel "
//					<< devFillLevels[count].id << devFillLevels[count].level <<  std::endl;
			++count;

		}

		/* Sort by increasing fill level */
		std::sort(fillLevels, fillLevels + numWorkers);

//		count = 0;
//		std::cout << "Dispatcher.cpp :" << "numItems " << req->size() <<  std::endl;
//		std::cout << "Dispatcher.cpp :" << "Push Item " << count++
//							<<  " to Device " << devFillLevels[0].id <<  std::endl;

		/* push first item to Worker with the lowest level first */
		pushItemToCPUWorker(req->nextItem(), fillLevels[0].id);

		/* The following statements conserve that devFillLevels stays sorted */
		for (unsigned i = 0; i < req->size()-1; ++i) {
			for (unsigned j = 0; j < numWorkers-1; ++j) {
				if(fillLevels[j] < fillLevels[j + 1]) {

//					std::cout << "Dispatcher.cpp :" << "0 Push Item " << count++
//							<<  " to Device " << devFillLevels[j].id <<  std::endl;

					pushItemToCPUWorker(req->nextItem(), fillLevels[j].id);
					++fillLevels[j].level;
					break;
				}
				if (j == numWorkers-2) {

//					std::cout << "Dispatcher.cpp :" << "1 Push Item " << count++
//							<<  " to Device " << devFillLevels[numWorkers-1].id <<  std::endl;

					pushItemToCPUWorker(req->nextItem(), fillLevels[numWorkers-1].id);
					++fillLevels[numWorkers-1].level;
				}
			}
		}
	}

}

int as::Dispatcher::createRequest(DOF* dofs, const unsigned& numDOFs,
		const int& gridId, const int& recId, const int& ligId,
		const Request::useMode_t& mode)
{
	/* check validity only in Debug mode */
	assert(gridValid(gridId));
	assert(protValid(recId));
	assert(protValid(ligId));

	int numItems = (numDOFs + _itemSize - 1) / _itemSize;
	WorkerItem* items[numItems];
	/* fetch items from buffer management */
	for (int i = 0; i < numItems; ++i) {
		items[i] = getItemFromBufMngt(); // Is thread safe!
		//DEBUG
		if (items[i] == NULL) {
			/* return valid items to buffer management */
			for (int j = 0; j < i; ++j) {
				assert(items[j] != NULL);
				returnItemToBufMngt(items[j]);
			}
			return -1;
		}
	}
	/* copy data to items */
	for (int i = 0; i < numItems-1; ++i) {
		WorkerItem* item = items[i];
		std::memcpy(item->DOFBuffer(), dofs + i*_itemSize, _itemSize*sizeof(DOF));
		item->setNumDOFs(_itemSize);
		item->setGlobGridId(gridId);
		item->setGlobRecId(recId);
		item->setGlobLigId(ligId);
	}
	int numLast = numDOFs - (numItems-1) * _itemSize;
	WorkerItem* item = items[numItems-1];
	std::memcpy(item->DOFBuffer(), dofs + (numItems-1)*_itemSize, numLast*sizeof(DOF));
	item->setNumDOFs(numLast);
	item->setGlobGridId(gridId);
	item->setGlobRecId(recId);
	item->setGlobLigId(ligId);

	/* assigned Id and push Items to request */
	int id = _id++; // atomically
	/* assert that request with this id is not yet contained */
	assert(!reqIdValid(id)); // request with this id is not yet contained;
	// insert a new id and get the reference to default created object

	_mutexMap.lock();
	Request* req = new Request();
	_reqMap[id] = req;
	_mutexMap.unlock();

	for (int i = 0; i < numItems; ++i) {
		req->addItem(items[i]);
	}
	req->setUseMode(mode);

	{
		std::lock_guard<std::mutex> guard(_mutexQueue);
		_reqQueue.push(req);
		_condVar.notify_one();
	}

	return id;
}

as::Dispatcher::pullErr_t as::Dispatcher::pullRequest(const int& id, EnGrad* const buffer)
{
	Request* req;
	{
		std::lock_guard<std::mutex> guard(_mutexMap);
		/* check if request has already been processed */
		if (!reqIdValid(id)) {
			return invalid;
		}
		/* no request is create because we check validity */
		req = _reqMap[id];
		if (!req->ready()) {
			return notReady;
		} else {
			/* invalidate id and go on */
			_reqMap.erase(id);
		}
	}
	unsigned size = req->size();
	unsigned shift = 0;
	for (unsigned i = 0; i < size; ++i) {
		WorkerItem* item = req->removeItem();
		unsigned itemSize = item->size();
		std::memcpy(buffer + shift, item->EnGradBuffer(), itemSize*sizeof(EnGrad));
		shift += itemSize;

		/* recyle item */
		item->reset();
		returnItemToBufMngt(item);
	}
	delete req;
	return ready;
}

/*
 ** @brief: this is done by the dispatcher itself
 */
inline as::Request* as::Dispatcher::removeRequest() {
	std::unique_lock<std::mutex> uniqueLock(_mutexQueue);
	while(_reqQueue.size() == 0) {
		if (_terminate.load() == true) {
			_sleeping = false;
//				std::cout << "Dispatcher.h :" << "Dispatcher goes to sleep" << std::endl;
			return NULL;
		}
		_sleeping = true;
		_condVar.wait(uniqueLock);
	}
//		std::cout << "Dispatcher.h :" << "Dispatcher wakes up" << std::endl;
	_sleeping = false;
	Request* req = _reqQueue.front();
	_reqQueue.pop();
	return req;
}


/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/

inline as::WorkerItem* as::Dispatcher::getItemFromBufMngt() {
		return _S_mngt.BufferMngt()->removeItem();
	}

inline void as::Dispatcher::returnItemToBufMngt(WorkerItem* const item) {
	_S_mngt.BufferMngt()->addItem(item);
}

inline bool as::Dispatcher::gridValid(const int& globalId) {
	return _S_mngt.DataMngt()->gridUnionValid(globalId);
}

inline bool as::Dispatcher::protValid(const int& globalId) {
	return _S_mngt.DataMngt()->proteinValid(globalId);
}

inline const std::set<unsigned>& as::Dispatcher::getCommonDeviceIDs_cref (const int& globalGridId, const int& globalProteinId0, const int& globalProteinId1) const {
	return _S_mngt.DataMngt()->commonDeviceIDs(globalGridId, globalProteinId0, globalProteinId1);
}

inline unsigned as::Dispatcher::deviceFillLevel(const unsigned deviceId) const {
	return _S_mngt.WorkerMngt()->deviceFillLevel(deviceId);
}

inline unsigned as::Dispatcher::hostFillLevel(const unsigned id) const {
	return _S_mngt.WorkerMngt()->hostFillLevel(id);
}

inline void as::Dispatcher::pushItemToGPUWorker(WorkerItem* const item, const unsigned& deviceId) {
	_S_mngt.WorkerMngt()-> pushItemToGPUWorker(item, deviceId);
}

inline void as::Dispatcher::pushItemToCPUWorker(WorkerItem* const item, const unsigned& id) {
	_S_mngt.WorkerMngt()-> pushItemToCPUWorker(item, id);
}

inline unsigned as::Dispatcher::numCPUWorker() const {
	return _S_mngt.WorkerMngt()->numCPUWorkers();
}
