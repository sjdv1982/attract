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
#include <cassert>
#include <mutex>

#include "as/BufferManagement.h"
#include "asUtils/macros.h"
#include "asUtils/Logger.h"

using namespace Log;


/* Constructor */
as::BufferManagement::BufferManagement(unsigned numItems, unsigned numDOFsPerItem) :
		_numItems(numItems),
		_contDOFBuffer(nullptr),
		_contEnGradBuffer(nullptr),
		_contItemBuffer(nullptr),
		_itemQueue()
{
	assert(numItems > 0);
	assert(numDOFsPerItem > 0);
	assert(numItems*numDOFsPerItem*(sizeof(DOF)+sizeof(EnGrad)) <
			1000*1000*10*(sizeof(DOF)+sizeof(EnGrad))); // < 540 MByte == 10 mio DOFs

	cudaError e = cudaMallocHost((void**)&_contDOFBuffer, numItems*numDOFsPerItem*sizeof(DOF));
	if (e == cudaSuccess) {
		CUDA_CHECK( cudaMallocHost((void**)&_contEnGradBuffer, numItems*numDOFsPerItem*sizeof(EnGrad)) );
		ASSERT(_contDOFBuffer != nullptr && _contEnGradBuffer != nullptr);
		cudaMallocHostused = true;
	} else {
		ASSERT(_contDOFBuffer == nullptr);
		_contDOFBuffer = new DOF[numItems*numDOFsPerItem];
		_contEnGradBuffer = new EnGrad[numItems*numDOFsPerItem];
		cudaMallocHostused = false;
	}

	_contItemBuffer = new WorkerItem[numItems];

	/* Cut large buffer into small pieces to create work items
	 * Store the Item in the queue container */
	for (unsigned i = 0; i < numItems; ++i) {
		as::WorkerItem* item = _contItemBuffer + i;
		item->setDOFBuffer(_contDOFBuffer + i*numDOFsPerItem);
		item->setEnGradBuffer(_contEnGradBuffer + i*numDOFsPerItem);
		_itemQueue.push(item);
	}
}

/* Destructor
 * Make sure the no items are currently processed
 * or have not yet been fetched by the client.
 * For that call the complete() method */
as::BufferManagement::~BufferManagement() {
	global_log->info() << "Shutdown Buffer-Management:" << std::endl;
	if (!complete()) {
		global_log->warning() << "Warning: Item buffers are getting freed "
				<< "while computations might still be in progress." << std::endl;
	} else {
		global_log->info() << std::setw(5) << " " << "Ok" << std::endl;
	}

	if (cudaMallocHostused == true) {
		CUDA_CHECK(cudaFreeHost(_contDOFBuffer));
		CUDA_CHECK(cudaFreeHost(_contEnGradBuffer));
	} else {
		delete[] _contDOFBuffer;
		delete[] _contEnGradBuffer;
	}
	delete[] _contItemBuffer;
}



/****************************
 * public member functions
 ****************************/
as::WorkerItem* as::BufferManagement::removeItem()
{
	std::lock_guard<std::mutex> guard(_mutex);
	if (_itemQueue.size() == 0) {
		return NULL;
	}

	WorkerItem* item = _itemQueue.front();
	_itemQueue.pop();
	return item;
}

void as::BufferManagement::addItem(WorkerItem* const item)
{
	std::lock_guard<std::mutex> guard(_mutex);
	assert(item->isReady() == false);
	_itemQueue.push(item);
}

/*
 ** @brief: Checks if all items are back to the buffer management.
 */
bool as::BufferManagement::complete()
{
	// maybe protection is unnecessary
	std::lock_guard<std::mutex> guard(_mutex);
	return _itemQueue.size() == _numItems;
}



/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


