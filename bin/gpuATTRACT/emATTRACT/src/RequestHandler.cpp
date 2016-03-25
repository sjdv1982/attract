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
#include <cassert>
#include <algorithm>

#include <AttractServer>
#include "nvToolsExt.h"

#include "RequestHandler.h"
#include "Chunk.h"
#include "SolverFactoryImpl.h"

using std::cerr;
using std::cout;
using std::endl;


void ema::RequestHandler::init(extServer& server, std::string const& solverName, std::vector<extDOF>& dofs) {
	_server = &server;
	std::unique_ptr<SolverFactory> factory(new SolverFactoryImpl);

	/* Fill the object array */
	for (unsigned i = 0; i < dofs.size(); ++i) {
		SharedSolver solver;
		solver = factory->createSolverByName(solverName);

		solver->setState(dofs[i]);
		_objects.emplace_hint(_objects.end(), i, solver);
	}

	/* set number of recieved objects */
	_numObjects = _objects.size();

	/*
	 * initialize chunk list based on the number of available structures
	 */


	/* shrink the number of chunks to fit the minimal chunkSize */
	while (_numObjects < _numChunks*_minChunkSize && _numChunks > 1) {
		 --_numChunks;
	}
	assert(_numChunks >= 1);

	/* calculate chunk sizes */
	unsigned base_numPerChunk = MIN(_numObjects,_numConcurrentObjects) / _numChunks;
	unsigned rest = MIN(_numObjects,_numConcurrentObjects) % _numChunks;

	unsigned chunkSizes[_numChunks];
	std::fill(chunkSizes, chunkSizes + _numChunks, base_numPerChunk);

	assert(rest < _numChunks);
	for (unsigned i = 0; i < rest; ++i) {
		++chunkSizes[i];
	}

	/* setup chunks and put them into chunk list */
	for (unsigned i = 0; i < _numChunks; ++i) {
		_chunkList.emplace_back();
	}

	unsigned count = 0;
	for (auto& chunk : _chunkList) {
		ObjMapIter mapIter = _objects.begin();
		for (unsigned i = 0; i < chunkSizes[count]; ++i, ++mapIter) {
			chunk.getContainer().push_back(std::move(*mapIter));
			_objects.erase(mapIter);
			assert(mapIter != _objects.end());
		}
		assert(count < _numChunks);
		++count;

	}
}

//#define H_IO

void ema::RequestHandler::run() {

	_collectedRequests.reserve(_chunkList.begin()->size());
	_collectedResults.reserve(_chunkList.begin()->size());

//	RingArray<int> reqIds;

	/* initial loop: start solvers and collect first requests and submit*/
	for (auto& chunk : _chunkList) {
		nvtxRangePushA("Processing");
		for (auto& obj : chunk.getContainer()) {
			SharedSolver& solver = obj.second;
			solver->start();
			_collectedRequests.push_back(Vector2extDOF(solver->getState()));
		}
		nvtxRangePop();

		nvtxRangePushA("Submit");
//		int reqId = ema::server_submit(*_server, _collectedRequests.data(), _collectedRequests.size(),
//				_serverOpt.gridId, _serverOpt.recId, _serverOpt.ligId, _serverOpt.useMode);
		int reqId = asClient::server_submit(*_server, _collectedRequests.data(), _collectedRequests.size(),
				_serverOpt.gridId, _serverOpt.recId, _serverOpt.ligId, _serverOpt.useMode);
		nvtxRangePop();
		chunk.setFetchSize(chunk.size());

		if (reqId == -1) {
			cerr << "Error: Submitting request." << std::endl;
			std::exit(EXIT_FAILURE);
		}
//		reqIds.push_back(reqId);
		chunk.setReqId(reqId);

		_collectedRequests.resize(0);
	}


	unsigned count = 1;
	while(_finishedObjects.size () < _numObjects && count < 100000) {

		/* load balancing */

		/* Adjust chunk sizes to balance the workload of each chunk.
		 * This happens in case that the global object list is empty and
		 * the chunks cannot be refilled by new initial configurations.
		 * Do it not each iteration */
		if (_objects.empty() && count%4 == 0) {
			double ratio = chunkSizeRatio(_chunkList);
			if(ratio > 1.5) {
				loadBalanceChunks(_chunkList);
			}
		}


			for (auto chunkListIter = _chunkList.begin(); chunkListIter != _chunkList.end(); ) {
				auto& chunk = *chunkListIter;
				_collectedRequests.resize(0);
				if (chunk.fetchSize() > 0) {
					_collectedResults.resize(chunk.fetchSize());


					/* Wait for requests */
					nvtxRangePushA("Waiting");
	//				unsigned count = ema::server_pull(*_server, chunk.reqId(), _collectedResults.data());
					unsigned count = asClient::server_pull(*_server, chunk.reqId(), _collectedResults.data());
					nvtxRangePop();

					if (count >= 10000) {
						cerr << "Error: pulling for Request." << std::endl;
						std::exit(EXIT_FAILURE);
					}

					/* Assigne results */
					chunk.setResults(_collectedResults);
				}

			/* Check if other chunks assigned results (after loadbalancing)*/
			chunk.checkLBconts();

			/* Process chunk and remove converged structures */
			nvtxRangePushA("Processing");
			auto iter = chunk.getContainer().begin();
			iter = chunk.getContainer().begin();

			while (iter != chunk.getContainer().end()) {
				SharedSolver& solver = iter->second;
				solver->step();

				/* test for convergence */
				if(solver->converged()) {

					/* destroy coroutine context by calling finalize */
					solver->finalize();
					/* move structure/object in finished object container */
					_finishedObjects.insert(move(*iter)); // no copy-construction
					chunk.getContainer().erase(iter++);

					/* move new structure/solver from object map if any left*/
					if (!_objects.empty()) {

						ObjMapIter objIter = _objects.begin();
						iter = chunk.getContainer().insert(iter, std::move(*objIter));
						_objects.erase(objIter);

						/* prepare new solver */
						SharedSolver& newSolver = iter->second;
						newSolver->start();
						/* collect new request */
						_collectedRequests.push_back(Vector2extDOF(newSolver->getState()));
						++iter;
					}
				} else {
					/* collect new request */
					_collectedRequests.push_back(Vector2extDOF(solver->getState()));
					++iter;
				}

			}
			nvtxRangePop();
			assert(iter == chunk.getContainer().end());

			chunk.setFetchSize(chunk.size());

			/* submit request */
			if (_collectedRequests.size() > 0) { // there is still something to submit
				nvtxRangePushA("Submit");
//				int reqId = ema::server_submit(*_server, _collectedRequests.data(), chunk.size(),
//						_serverOpt.gridId, _serverOpt.recId, _serverOpt.ligId, _serverOpt.useMode);
				int reqId = asClient::server_submit(*_server, _collectedRequests.data(), chunk.size(),
						_serverOpt.gridId, _serverOpt.recId, _serverOpt.ligId, _serverOpt.useMode);
				nvtxRangePop();

				if (reqId == -1) {
					cerr << "Error: Submitting request." << std::endl;
					std::exit(EXIT_FAILURE);
				}
				chunk.setReqId(reqId);
			}

			++chunkListIter;
			// do not remove any chunks since objects might still reside in LBconts waiting for results from other chunks

		} // for each chunk

		++count;
	} // while
	assert(_finishedObjects.size() == _numObjects);

}

std::vector<ema::extDOF> ema::RequestHandler::getResultStates() {
	std::vector<extDOF> stateVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		stateVec[i] = Vector2extDOF(_finishedObjects[i]->getState());
	}
	return stateVec;
}
std::vector<ema::extEnGrad> ema::RequestHandler::getResultEnGrads() {
	std::vector<extEnGrad> enGradVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		enGradVec[i] = ObjGrad2extEnGrad(_finishedObjects[i]->getObjective());
	}
	return enGradVec;
}
std::vector<std::unique_ptr<ema::Statistic>> ema::RequestHandler::getStatistics() {
	std::vector<std::unique_ptr<Statistic>> statisticVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		statisticVec[i] = _finishedObjects[i]->getStats();
	}
	return statisticVec;

}
