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

#ifndef SHARED_H_
#define SHARED_H_

#include <cassert>
#include <mutex>
#include <memory>

namespace as {

/*
 ** @brief: This class is for trivial and thread safe
 ** reference counting of a shared object.
 ** The counting itself is performed externally.
 ** The shared object is destroyed
 ** 	either the count decreases to zero
 ** 	or the shared variable is destroyed.
 */
template <class T>
class Shared {
public:
	/* Constructor */
	Shared(T* object) : _object(object), _count(1), _valid(true) {}
	Shared() : _object(NULL), _count(0), _valid(false) {}

	/* Destructor */
	~Shared() {}

	/***************
	* G E T T E R
	***************/
	T* get() const {
		assert(_valid == true);
		return _object.get();
	}

	bool isValid() const {
		return _valid;
	}

	unsigned useCount() const {
		return _count;
	}

	unsigned real_use_count() const {
		return _object.use_count();
	}

	/***************
	* S E T T E R
	***************/
	void set(T* object) {
		assert(_valid == false);
		assert(_count == 0);
		std::shared_ptr<T> ptr(object);
		_object = ptr;
		++_count;
		_valid = true;
	}

	/****************************
	 * public member functions
	 ****************************/
	void incrCount() {
		assert(_valid == true);
		assert(_count > 0);
		++_count;
	}

	void decrCount() {
		assert(_count > 0);
		--_count;
		if (_count == 0) {
			assert(_object.unique());
			_valid = false;
			_object.reset();
		}
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

	std::shared_ptr<T> _object;		/** shared pointer to contained object */
	unsigned _count;	/** number of referring (this object using) objects */
	bool _valid;
};


} // namespace AServ

#endif /* SHARED_H_ */
