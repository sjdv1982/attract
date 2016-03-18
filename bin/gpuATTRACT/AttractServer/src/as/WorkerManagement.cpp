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

#include "as/WorkerManagement.h"
#include "asUtils/Logger.h"

using namespace Log;

/* Constructor */
as::WorkerManagement::WorkerManagement(ServerManagement& S_mngt, unsigned numDOFsPerItem, unsigned GPU_atomBufferSize, unsigned CPU_atomBufferSize) :
		_GPU_AtomBufferSize(GPU_atomBufferSize),
		_CPU_AtomBufferSize(CPU_atomBufferSize),
		_dofBufferSize(numDOFsPerItem),
		_S_mngt(S_mngt) {}

/* Destructor */
as::WorkerManagement::~WorkerManagement() {
	global_log->info() << "Shutdown Worker-Management:" << std::endl;
	for (int i :_GPUWorkerIds) {
		global_log->info() << std::setw(5) << " " << "Signal to terminate GPUWorker " << i << std::endl;;
		GPUWorker* worker = _GPUWorkers[i];
		worker->signalTerminate();
		worker->join();
		delete worker;
	}

	for (auto worker : _CPUWorkers) {
		global_log->info() << std::setw(5) << " " << "Signal to terminate CPUWorker " << worker->id() << std::endl;

		worker->signalTerminate();
		worker->join();
		delete worker;
	}

	global_log->info() << std::setw(5) << " " << "Ok" << std::endl;
}



/****************************
 * public member functions
 ****************************/

unsigned as::WorkerManagement::numCPUWorkers() {
	std::lock_guard<std::mutex> guard(_m_CPUWorker);
	return _CPUWorkers.size();
}

void as::WorkerManagement::addGPUWorker(int deviceId) {
	/* assert that id is not occupied */
	assert(_GPUWorkerIds.find(deviceId) == _GPUWorkerIds.end());
	_GPUWorkerIds.insert(deviceId);

	GPUWorker* worker = new GPUWorker(_S_mngt, deviceId, _GPU_AtomBufferSize, _dofBufferSize);
	_GPUWorkers.placeAtLoc(worker, deviceId);
	worker->start();
}

void as::WorkerManagement::addCPUWorker() {
	std::lock_guard<std::mutex> guard(_m_CPUWorker);
	CPUWorker* worker = new CPUWorker(_S_mngt, _CPU_AtomBufferSize);
	_CPUWorkers.push_back(worker);
	worker->set(_CPUWorkers.size()-1);
	worker->start();
}

void as::WorkerManagement::removeGPUWorker(int deviceId) {
	/* assert that id is not occupied */
	assert(_GPUWorkerIds.find(deviceId) != _GPUWorkerIds.end());
	GPUWorker* worker = _GPUWorkers[deviceId];
	worker->signalTerminate();
	worker->join();
	_GPUWorkerIds.erase(deviceId);
	delete worker;
}

void as::WorkerManagement::removeCPUWorker() {
	/* remove last CPU worker */
	std::lock_guard<std::mutex> guard(_m_CPUWorker);
	CPUWorker* worker = *(_CPUWorkers.end()-1);
	_CPUWorkers.pop_back();
	worker->signalTerminate();
	worker->join();
	delete worker;
}


/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/



