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

#ifndef BUFFERMANAGEMENT_H_
#define BUFFERMANAGEMENT_H_

#include <list>
#include <mutex>
#include <condition_variable>
#include <queue>

#include <as/WorkItem.h>

namespace as {


class BufferManagement {
public:
	/* Constructor */
	BufferManagement(unsigned numItems, unsigned numDOFsPerItem);

	/* Destructor */
	~BufferManagement();

	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	WorkerItem* removeItem();

	void addItem(WorkerItem* const item);

	/*
	 ** @brief: Indicates if all items are back
	 ** in the queue. The object can now safely be deleted.
	 */
	bool complete();

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

	unsigned _numItems;

	DOF* _contDOFBuffer;		/** Page-Locked allocation */
	EnGrad* _contEnGradBuffer;	/** Page-Locked allocation */
	WorkerItem* _contItemBuffer;/** Ptr. to all allocated WorkerItems */

	bool cudaMallocHostused;

	std::queue<WorkerItem*, std::list<WorkerItem*> > _itemQueue;

	std::mutex _mutex;
};

} // namespace

#endif /* BUFFERMANAGEMENT_H_ */
