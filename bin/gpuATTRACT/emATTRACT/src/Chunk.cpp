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
#include "Chunk.h"

void ema::Chunk::setResults(const std::vector<extEnGrad>& results) {
	auto iter = results.begin();
	/* fill main list */
	for (auto& pair : _cont ) {
		SharedSolver& solver = pair.second;
		solver->setObjective(*iter);
		++iter;
	}

	/* fill LB lists, namely objectes contained by other chunks */
	for (auto LBcontPtr : _LBcontRefs) {
		for (auto& pair : LBcontPtr->cont) {
			SharedSolver& solver = pair.second;
			solver->setObjective(*iter);
			++iter;
		}
		LBcontPtr->complete = true;
	}
//	std::cerr << "setResults: " << _cont.size() << "     ";
//	for (auto LBcontPtr : _LBcontRefs) {
//		std::cerr << LBcontPtr->cont.size() << " ";
//	}
//	std::cerr << std::endl;
	assert(iter == results.end());

	/* remove LBcontRefs */
	_LBcontRefs.clear();
}

void ema::Chunk::checkLBconts() {
	for (auto LBcontIter = _LBconts.begin();  LBcontIter != _LBconts.end();) {

		if (LBcontIter->complete) {
			for (auto& ShSolver : LBcontIter->cont) {
				_cont.push_back(std::move(ShSolver));
			}
			/* destroy container */
			_LBconts.erase(LBcontIter++);
		} else {
			++LBcontIter;
		}

	}
}

unsigned ema::Chunk::overAllSize() {
	unsigned count = 0;
	count += _cont.size();
//	std::cerr << "overAllSize: " << _cont.size() << "     ";
	for (auto& LBCont : _LBconts) {
		count += LBCont.cont.size();
//		std::cerr << LBCont.cont.size() << " ";
	}
//	std::cerr << std::endl;
	return count;
}

//#define H_IO

void ema::loadBalanceChunks (std::list<Chunk>& chunkList) {
	using std::cerr;
	using std::endl;

	assert(chunkList.size() > 1);
	unsigned numObjects = 0;
	for (auto& chunk : chunkList) {
		assert(chunk.size() == chunk.overAllSize());
		numObjects += chunk.size();
	}

	unsigned base_numPerChunk = numObjects / chunkList.size();
	unsigned rest = numObjects % chunkList.size();

	unsigned chunkSizes[chunkList.size()]; // required chunk sizes
	std::fill(chunkSizes, chunkSizes + chunkList.size(), base_numPerChunk);

	assert(rest < chunkList.size());
	for (unsigned i = 0; i < rest; ++i) {
		++chunkSizes[i];
	}

	/* if one of the chunk sizes is below a threshold we erase that chunk and distribute its contents to the remaining chunks */
	/* does not work, since we still need to wait for the results of the erased chunk */
//	unsigned i = 0;
//	for (auto chunk_iter = chunkList.begin(); chunk_iter != chunkList.end(); ++i) {
//		if (chunkSizes[i] < _minChunkSize) {
//
//			auto chunk_iter_to = (chunk_iter != chunkList.begin() ? chunkList.begin() : --chunkList.end());
//
//		} else {
//			++chunk_iter;
//		}
//	}

	/* generate a temp vector to sort the chunks according to their size in descending order */
	std::vector<Chunk*> chunkVec(chunkList.size());

	int tmp = 0;
	for (auto& chunk : chunkList) {
		chunkVec[tmp++] = &chunk;
	}

#ifdef H_IO
	cerr << "Chunk sizes before sorting" << endl;
	for (auto& chunk : chunkVec) {
		cerr << chunk->size() << " ";
	}
	cerr << endl;
#endif
	/* sort in descending order */
	std::sort(chunkVec.begin(), chunkVec.end(), [](const Chunk* lhs, const Chunk* rhs) -> bool { return lhs->size() > rhs->size();});

#ifdef H_IO
	cerr << "Chunk sizes after sorting" << endl;
	for (auto& chunk : chunkVec) {
		cerr << chunk->size() << " ";
	}
	cerr << endl;
#endif


	unsigned i = 0;
	for (auto chunk_iter_from = chunkVec.begin(); chunk_iter_from != (--chunkVec.end()); ++chunk_iter_from, ++i) {
		int surplus = (*chunk_iter_from)->overAllSize() - chunkSizes[i];
//		cerr << " chunkSizes[" << i << "]="<< chunkSizes[i] << " (*chunk_iter_from)->size()=" << (*chunk_iter_from)->size()
//					<< " surplus=" << surplus << endl;
		assert(surplus >= 0);
		/* distribute the elements to chunks that contain less elements to the right */
		auto chunk_iter_to = chunk_iter_from;
		for (int j = i+1; surplus > 0; j++) {
			++chunk_iter_to;
			assert(chunk_iter_to != chunkVec.end());
			const int deficit = chunkSizes[j] - (*chunk_iter_to)->overAllSize();
			if (deficit <= 0) continue;

			/* since we are now moving elements to another chunk we need to set its fetchSize() to its current size.
			 * otherwise to less elements are pulled from the server
			 */
			(*chunk_iter_from)->setFetchSize((*chunk_iter_from)->size());
//			(*chunk_iter_to)->setFetchSize((*chunk_iter_to)->size());

#ifdef H_IO
			cerr << "deficit=" << deficit << " chunkSizes[" << j << "]="<< chunkSizes[j] << " (*chunk_iter_to)->size()=" << (*chunk_iter_to)->size()
					<< " surplus=" << surplus << endl;
			cerr << "moving " << std::min(deficit,surplus) << " objects from chunk " << i << " to " << j << endl;
#endif
			int min_max = std::min(deficit,surplus);
			auto& LBcont = (*chunk_iter_to)->createLBCont(min_max);
			for (int k = 0; k < min_max; ++k) {
				LBcont.cont[LBcont.cont.size()-1 - k] = std::move((*chunk_iter_from)->getContainer().back()); // move-assignment
				(*chunk_iter_from)->getContainer().pop_back(); // remove transfered item; note object of coro (ptr) gets not destroyed
				// since move-assignment sets it to null.
			}
			(*chunk_iter_from)->acceptLBCont(LBcont);
			surplus -= min_max;
		} // for
	} // for
}

double ema::chunkSizeRatio (std::list<Chunk> const& chunkList) {

	auto fnc = [] (Chunk const& lhs, Chunk const& rhs) -> bool { return lhs.size() < rhs.size();} ;
	auto minmax = std::minmax_element(chunkList.begin(), chunkList.end(), fnc);
	auto min = minmax.first->size();
	auto max = minmax.second->size();
	if ((max <= 2 && min == 1) || max <= 1)  {
		return 0.0;
	}

	return static_cast<double>(max) / static_cast<double>(min);


}

