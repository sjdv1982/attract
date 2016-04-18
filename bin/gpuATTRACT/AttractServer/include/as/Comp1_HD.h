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

#ifndef COMP1_HD_H_
#define COMP1_HD_H_

#include "as/cudaHostAllocation.h"

namespace as {

template <class T, class ALLOC>
class Comp1_HD {
public:
	/* Constructor */
	Comp1_HD (unsigned numEl) :
	_numEl(numEl),
	_h_data(NULL),	_d_data(NULL),
	_cpySize(_numEl*sizeof(T)),
	_deviceAlloc(false), _hostAlloc(false){ }

	/* Destructor */
	~Comp1_HD() {
		freeDevice();
		freeHost();
	}


	/***************
	* G E T T E R
	***************/
	inline unsigned size() const {
		return _numEl;
	}

	inline T* h_data() const {
		return _h_data;
	}

	inline T* d_data() const {
		return _d_data;
	}

	/***************
	* S E T T E R
	***************/

	inline void set_h_data(T* data) {
		_h_data = data;
	}

	inline void set_d_data(T* data) {
		_d_data = data;
	}

	/****************************
	 * public member functions
	 ****************************/
	void initHost() {
		if (!_hostAlloc && ALLOC::HostAlloc) {
			ALLOC::malloc(_h_data, _numEl*sizeof(T));
			memset(_h_data, 0, _numEl*sizeof(T));
			_hostAlloc = true;
		}
	}

	void initDevice() {
		if (!_deviceAlloc && ALLOC::DeviceAlloc) {
			cudaVerify(cudaMalloc((void**)&_d_data, _numEl*sizeof(T) ));
			cudaVerify(cudaMemset(_d_data, 0, _numEl*sizeof(T) ));
			_deviceAlloc = true;
		}
	}

	void freeDevice() {
		if (_deviceAlloc) {
			cudaFree(_d_data);
			_deviceAlloc = false;
		}
	}

	void freeHost() {
		if (_hostAlloc) {
			ALLOC::free(_h_data);
			_hostAlloc = false;
		}
	}

	void resetHostData() {
		if (_hostAlloc) {
			memset(_h_data, 0, _numEl*sizeof(T));
		}
	}

	void resetDeviceData(const cudaStream_t &stream = 0) {
		if (_deviceAlloc) {
			cudaVerify(cudaMemsetAsync(_d_data, 0, _numEl*sizeof(T) , stream));
		}
	}

	void resetData(const cudaStream_t &stream = 0) {
		resetHostData();
		resetDeviceData(stream);
	}

	inline __host__ void cpyH2D(const cudaStream_t &stream = 0) {
		cudaVerify(cudaMemcpyAsync(_d_data, _h_data, _cpySize, cudaMemcpyHostToDevice, stream));
	}

	inline __host__ void cpyD2H(const cudaStream_t &stream = 0) {
		cudaVerify(cudaMemcpyAsync(_h_data, _d_data, _cpySize, cudaMemcpyDeviceToHost, stream));
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

	unsigned _numEl;		/** number of elements */

	T *_h_data; 		/** Host data, x-component */
	T *_d_data; 		/** Device data, x-component */

	unsigned _cpySize;		/** size of memory transactions in bytes */
	bool _deviceAlloc;	/** memory at device/host already allocated? */
	bool _hostAlloc;
};

}

#endif /* COMP1_HD_H_ */
