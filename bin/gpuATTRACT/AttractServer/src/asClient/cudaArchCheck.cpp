/*
 * cudaArchCheck.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: uwe
 */

#include "cuda_runtime.h"
#include <exception>
#include <string>
#include <sstream>
#include "asUtils/macros.h"

#include <iostream>

constexpr int minMajorCC = 3;
constexpr int minMinorCC = 0;

class CUDAComputeCapabilityException : public std::exception {
public:
	CUDAComputeCapabilityException(int major, int minor, int deviceId) :
		_majorCC(major),
		_minorCC(minor),
		_deviceId(deviceId)
	{}

	virtual const char* what() const throw()
	{
		std::stringstream msg;
		msg << "CUDA Compute Capability ("
				<< _majorCC << "." << _minorCC << ") "
				<< "of device (ID) " << _deviceId << " is insufficient. At least "
				<< minMajorCC << "." << minMinorCC << " is required.";
		return msg.str().c_str();
	}

private:
	int _majorCC;
	int _minorCC;
	int _deviceId;

};

namespace asClient {

void checkComputeCapability() {
	int nDevices = 0;
	CUDA_CHECK(cudaGetDeviceCount(&nDevices));


	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		CUDA_CHECK(cudaGetDeviceProperties(&prop, i));

		int majorCC = prop.major;
		int minorCC = prop.minor;

		int CC = 10*majorCC + minorCC;
		int minCC = 10*minMajorCC + minMinorCC;

		if (minCC > CC) {
			throw CUDAComputeCapabilityException(majorCC, minorCC, i);
		}
	}
}

}


