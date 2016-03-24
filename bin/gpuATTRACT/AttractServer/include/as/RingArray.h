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

#ifndef RINGARRAY_H_
#define RINGARRAY_H_

namespace as {

template <class T>
class RingArray {
	public:
	/* Constructor */
	RingArray(unsigned size) : _size(size), _pointer(-1)
	{
		_ring = new T[_size];
	}

	/* Destructor */
	~RingArray() {
		delete[] _ring;
	}

	/* operator */
//	inline T& operator[] (const unsigned& stage) {
//		return _ring[calcIdx(stage)];
//	}

	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/

	inline void push(const T& obj) {
		rotate();
		_ring[calcIdx(0)] = obj;
	}

	inline T get(const unsigned& stage) const {
		return _ring[calcIdx(stage)];
	}


	inline void rotate() {
		_pointer = (_pointer + 1) % _size;
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

	inline unsigned calcIdx (const unsigned& stage) const {
		return ((int)_size -_pointer + (int)stage) % _size;
	}

	/****************************
	 * private member variables
	 ****************************/

	unsigned _size;
	T* _ring;
	int _pointer;

};

}  // namespace AServ


#endif /* RINGARRAY_H_ */
