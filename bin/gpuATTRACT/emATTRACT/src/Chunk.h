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

#ifndef CHUNK_H_
#define CHUNK_H_

#include <memory>
#include <list>
#include "SolverBase.h"
 
namespace ema {

class Chunk {
public:

	int reqId() const {
		return _reqId;
	}

	void setReqId(int value) {
		_reqId = value;
	}

	unsigned fetchSize() const {
		return _fetchSize;
	}

	void setFetchSize(unsigned value) {
		_fetchSize = value;
	}

	unsigned size() const {
		return _cont.size();
	}

	using SharedSolver = std::shared_ptr<SolverBase>;
	using ContainerType = std::list<std::pair<unsigned, SharedSolver>>;
	ContainerType& getContainer() {
		return _cont;
	}

	using iterator = ContainerType::iterator;

	/*
	 * Container used to store objects from other chunks after a load balancing step.
	 * It is filled by other chucks that hold a reference/pointer to it thereby
	 * setting the complete flag to true.
	 * When *this chunk is processed, it checks if the complete flag is set. In this
	 * case, the elements of LBCont are moved to the main container (_cont).
	 *
	 */
	class LBCont {
	public:

		LBCont(unsigned numEl) : cont(numEl), complete(false) {}
		LBCont() : LBCont(0) {}

		std::vector<std::pair<unsigned, SharedSolver>> cont;
		bool complete;
	};


	/* is called by other chunks to get a LBCont they want to process.
	 * process means: the they assign the results and set the complete flag */
	LBCont& createLBCont(unsigned numEl) {
		_LBconts.push_back(numEl);
		return _LBconts.back();
	}

	/*
	 * accept a LBCont
	 */
	void acceptLBCont(LBCont& lbcont) {
		_LBcontRefs.push_front(&lbcont); // append at front to preseve original ordering of main list (_cont)
	}

	void setResults (const std::vector<extEnGrad>& results);

	/*
	 * checks if a LBCont in _LBconts is already processed. If so, the elements of the
	 * respective LBCont is appended to the main list (_cont) and is destroyed
	 * afterwards.
	 * Note: This function may not be called unil all requests have been fetched.
	 * After calling this function solver steps may be performed.
	 */
	void checkLBconts();

	/*
	 * Returns the expected size after a load balancing step. This includes the size of
	 * the main list + the number of objects in the LB containers in _LBconts
	 */
	unsigned overAllSize();

private:
	int _reqId = -1; 		/** id for the server request */
	unsigned _fetchSize; 	/** number of requests to fetch */
	ContainerType _cont; 	// we use a list to be able to erase arbitrary elements


	std::list<LBCont> _LBconts; 	/** container of LBcontainers that needs to get processed by other chunks */
	std::list<LBCont*> _LBcontRefs; /** container of LBcontainers that need to get processed by this chunks.
	 	 	 	 	 	 	 	 	 ** processed means: results get assigned to these containers */

};

//#define H_IO

/*
 * Performs a load balancing step of the Chunks. This function may not be called in subsequent iterations.
 */
void loadBalanceChunks (std::list<Chunk>& chunkList);


/*
 * Returns the ratio between the largest and the smallest chunk size (> 1.0)
 */
double chunkSizeRatio (std::list<Chunk> const& chunkList);

} // namespace


#endif /* CHUNK_H_ */
