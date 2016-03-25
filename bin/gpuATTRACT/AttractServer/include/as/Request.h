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

#ifndef REQUEST_H_
#define REQUEST_H_

#include <list>
#include <cassert>

#include "as/WorkItem.h"

namespace as {


class Request {
public:
	/* Types */
	enum  useMode_t{
		GPU,
		CPU,
		GPUCPU,	// currently not supported
		unspecified
	};

	/* Constructor */
	Request() : _useMode(unspecified) {}

	/* Destructor */
	~Request() {}



	/***************
	* G E T T E R
	***************/

	inline unsigned size() const {
		return _workItems.size();
	}

	inline useMode_t useMode() const {
		return _useMode;
	}

	inline unsigned globGridId() const {
		assert(_workItems.size() > 0);
		return _workItems.front()->globGridId();
	}

	inline unsigned globRecId() const {
		assert(_workItems.size() > 0);
		return _workItems.front()->globRecId();
	}

	inline unsigned globLigId() const {
		assert(_workItems.size() > 0);
		return _workItems.front()->globLigId();
	}

	/***************
	* S E T T E R
	***************/
	inline void setUseMode(const useMode_t& mode) {
		_useMode = mode;
	}

	/****************************
	 * public member functions
	 ****************************/
	inline bool ready() const {
		for(const WorkerItem* item : _workItems) {
			if (!item->isReady()) {
				return false;
			}
		}
		return true;
	}

	inline void addItem(WorkerItem* item) {
		_workItems.push_back(item);
	}

	inline WorkerItem* nextItem() {
		assert(_workItems.size() > 0);
		assert(_it != _workItems.end());
		WorkerItem* item = *_it;
		++_it;
		return item;
	}

	/*
	 ** @brief: Call before getting Items by calling nextItem().
	 */
	inline void resetIterator() {
		_it = _workItems.begin();
	}

	/*
	 ** @brief: should only be called after the result was delivered to the client.
	 ** One uses it to recycle the items.
	 */
	inline WorkerItem* removeItem() {
		assert(_workItems.size() > 0);
		WorkerItem* item = _workItems.front();
		_workItems.pop_front();
		return item;
	}

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
	useMode_t _useMode;

	std::list<WorkerItem*>  _workItems;
	std::list<WorkerItem*>::iterator _it;
};


}  // namespace AServ

#endif /* REQUEST_H_ */
