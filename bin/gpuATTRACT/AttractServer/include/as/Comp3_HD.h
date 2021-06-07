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

#ifndef COMP3_HD_H_
#define COMP3_HD_H_

#include "as/cudaHostAllocation.h"

namespace as {


template <class T, class ALLOC>
class Comp3_HD {
public:
	/* Constructor */
	Comp3_HD (unsigned numEl) :
	_numEl(numEl),
	_numAlloc(((_numEl)/WARP_SIZE + 1)*WARP_SIZE),
	_h_conf(NULL), _h_x(NULL),	_h_y(NULL),	_h_z(NULL),
	_d_conf(NULL), _d_x(NULL),	_d_y(NULL),	_d_z(NULL),
	_cpySize(_numAlloc*3*sizeof(T)),
	_cpySize2(_numAlloc*sizeof(unsigned short)),
	_deviceAlloc(false), _hostAlloc(false) { }

	/* Destructor */
	~Comp3_HD() {
		freeDevice();
		freeHost();
	}


	/***************
	* G E T T E R
	***************/
	inline unsigned size() const {
		return _numEl;
	}
	inline unsigned short* h_conf () const {
		return _h_conf;
	}
	inline T* h_x () const {
		return _h_x;
	}
	inline T* h_y () const {
		return _h_y;
	}
	inline T* h_z () const {
		return _h_z;
	}

	inline unsigned short* d_conf () const {
		return _d_conf;
	}
	inline T* d_x () const {
		return _d_x;
	}
	inline T* d_y () const {
		return _d_y;
	}
	inline T* d_z () const {
		return _d_z;
	}

	/***************
	* S E T T E R
	***************/
	inline void set_h_conf (unsigned short *h) {
		 _h_conf = h;
	}
	inline void set_h_x (T* h) {
		 _h_x = h;
	}
	inline void set_h_y (T* h) {
		 _h_y = h;
	}
	inline void set_h_z (T* h) {
		 _h_z = h;
	}

	inline void set_d_conf (unsigned short *h) {
		 _d_conf = h;
	}
	inline void set_d_x (T* h) {
		 _d_x = h;
	}
	inline void set_d_y (T* h) {
		 _d_y = h;
	}
	inline void set_d_z (T* h) {
		 _d_z = h;
	}


	/****************************
	 * public member functions
	 ****************************/
	void initHost() {
		if (!_hostAlloc && ALLOC::HostAlloc) {
			ALLOC::malloc(_h_conf, _numAlloc*sizeof(unsigned short) );
			memset(_h_conf, 0, _numAlloc*sizeof(unsigned short));
			ALLOC::malloc(_h_x, _numAlloc*3*sizeof(T));
			memset(_h_x, 0, _numAlloc*3*sizeof(T));
			_h_y =  _h_x + _numAlloc;
			_h_z =  _h_x + 2*_numAlloc;
			_hostAlloc = true;
		}
	}

	void initDevice() {
		if (!_deviceAlloc && ALLOC::DeviceAlloc) {
			cudaVerify(cudaMalloc((void**)&_d_conf, _numAlloc*sizeof(unsigned short) ));
			cudaVerify(cudaMemset(_d_conf, 0, _numAlloc*sizeof(unsigned short) ));
			cudaVerify(cudaMalloc((void**)&_d_x, _numAlloc*3*sizeof(T) ));
			cudaVerify(cudaMemset(_d_x, 0, _numAlloc*3*sizeof(T) ) );
			_d_y =  _d_x + _numAlloc;
			_d_z =  _d_y + _numAlloc;
			_deviceAlloc = true;
		}
	}

	void freeDevice() {
		if (_deviceAlloc) {
			cudaFree(_d_conf);
			cudaFree(_d_x);
			_deviceAlloc = false;
		}
	}

	void freeHost() {
		if (_hostAlloc) {
			ALLOC::free(_h_conf);			
			ALLOC::free(_h_x);
			_hostAlloc = false;
		}
	}

	void resetHostData() {
		if (_hostAlloc) {
			memset(_h_conf, 0, _numAlloc*sizeof(unsigned short));
			memset(_h_x, 0, _numAlloc*3*sizeof(T));
		}
	}

	void resetDeviceData(const cudaStream_t &stream = 0) {
		if (_deviceAlloc) {
			cudaVerify(cudaMemsetAsync(_d_conf, 0, _numAlloc*sizeof(unsigned short) , stream));
			cudaVerify(cudaMemsetAsync(_d_x, 0, _numAlloc*3*sizeof(T) , stream));
		}
	}

	void resetData(const cudaStream_t &stream = 0) {
		resetHostData();
		resetDeviceData(stream);
	}
	
	inline __host__ void cpyH2D(const cudaStream_t &stream = 0) {
		cudaVerify(cudaMemcpyAsync(_d_conf, _h_conf, _cpySize2, cudaMemcpyHostToDevice, stream));
		cudaVerify(cudaMemcpyAsync(_d_x, _h_x, _cpySize, cudaMemcpyHostToDevice, stream));
	}

	inline __host__ void cpyD2H(const cudaStream_t &stream = 0) {
		cudaVerify(cudaMemcpyAsync(_h_conf, _d_conf, _cpySize2, cudaMemcpyDeviceToHost, stream));
		cudaVerify(cudaMemcpyAsync(_h_x, _d_x, _cpySize, cudaMemcpyDeviceToHost, stream));
	}

	__host__ void printEl(uint numEl, uint obj) {
		int colWidth = 15; // column WIDTH for output
		int precisionSetting = std::cerr.precision();
		std::ios::fmtflags flagSettings = std::cerr.flags();


		std::cerr.setf(std::ios::fixed | std::ios::showpoint | std::ios::showpos);
		std::cerr.precision(5);

		std::cerr << "VectorArray.h:\n";
		std::cerr << "Object " << obj << std::endl;
		std::cerr << std::left << std::setw(5) << "num" << ": "
						<< std::setw(colWidth) << "conf" << " "
						<< std::setw(colWidth) << "x" << " "
						<< std::setw(colWidth) << "y" << " "
						<< std::setw(colWidth) << "z" << std::endl;
		for (unsigned int i = 0; i < numEl; i++) {
			uint idx = i;
			std::cerr << std::left << std::setw(5) << i << ": "
					<< std::setw(colWidth) << _h_conf[idx] << " "
					<< std::setw(colWidth) << _h_x[idx] << " "
					<< std::setw(colWidth) << _h_y[idx] << " "
					<< std::setw(colWidth) << _h_z[idx] << std::endl;
		}
		std::cerr.precision(precisionSetting);
		std::cerr.flags(flagSettings);
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
	unsigned _numAlloc;		/** number of allocated elements per component (x,y,z,w,...)*/

	unsigned short *_h_conf;   /** Host data, conformers */
	T *_h_x; 		/** Host data, x-component */
	T *_h_y; 		/** Host data, y-component */
	T *_h_z; 		/** Host data, z-component */

	unsigned short *_d_conf;   /** Device data, conformers */
	T *_d_x; 		/** Device data, x-component */
	T *_d_y; 		/** Device data, y-component */
	T *_d_z; 		/** Device data, z-component */

	unsigned _cpySize;		/** size of memory transactions in bytes */
	unsigned _cpySize2;	/** size of conf memory transactions in bytes */
	bool _deviceAlloc;	/** memory at device/host already allocated? */
	bool _hostAlloc;
};

} //namespace

#endif /* COMP3_HD_H_ */

