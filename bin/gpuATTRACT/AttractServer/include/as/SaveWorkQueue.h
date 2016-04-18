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

#ifndef SAVEWORKQUEUE_H_
#define SAVEWORKQUEUE_H_

#include <list>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace as {


template <class T>
class SafeWorkQueue {
public:
	/* Constructor */
	SafeWorkQueue() : _sleeping(false), _terminate(false) {}


	/* Destructor */

	/***************
	* G E T T E R
	***************/
	
	inline bool sleeps() const {
		return _sleeping;
	}

	inline bool terminates() const {
		return _terminate;
	}

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	inline T removeCondBlocked () {
		std::unique_lock < std::mutex > uniqueLock(_mutex);
		while( (_queue.size() == 0) ) {
			if (_terminate.load() == true) {
				_sleeping = false;
//				std::cout << "SaveWorkQueue.h :" << "Worker goes to terminate" << std::endl;
				return T();
			}
//			std::cout << "SaveWorkQueue.h :" << "Worker goes to sleep" << std::endl;
			_sleeping = true;
			_condVar.wait(uniqueLock);
		}
//		std::cout << "SaveWorkQueue.h :" << "Worker wakes up" << std::endl;

		T item = _queue.front();
		_queue.pop_front();

		return item;
	}

	inline T removeUnblocked () {
		std::lock_guard<std::mutex> guard(_mutex);
		if (_queue.size() == 0) {
			return T();
		}
		T item = _queue.front();
		_queue.pop_front();

		return item;
	}

	inline void add(T item) {
		std::lock_guard<std::mutex> guard(_mutex);
		//DEBUG
//		std::cout << "SafeQr bef add2qu check validity: addr " << item  << " "<< item->DOFBuffer()[0] << std::endl;
		_queue.push_back(item);
		_condVar.notify_one();
	}

	inline void signalTerminate() {
		_terminate = true;
		_condVar.notify_one();
	}

	inline unsigned size() const {
		return _queue.size();
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

	std::list<T> _queue;
	std::mutex _mutex;
	std::condition_variable _condVar;

	std::atomic<bool> _sleeping;
	std::atomic<bool> _terminate;
};

} // namespace

#endif /* SAVEWORKQUEUE_H_ */
