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

#ifndef HELPER_H_
#define HELPER_H_

#include "asUtils/asUtilsTypes.h"
#include "asUtils/macros.h"

namespace asUtils {

inline double sizeConvert(unsigned long long sizeInByte, sizeType_t type = byte) {
	return double(sizeInByte) / type;
}

inline __host__ __device__ bool isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

inline __host__ __device__ unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

inline __host__ __device__ void getNumBlocksAndThreads(const int& n, const int& maxBlocks, const int& maxThreads,
		int &blocks, int &threads)
{
	threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
	blocks = (n + (threads * 2 - 1)) / (threads * 2);
	blocks = MIN(maxBlocks, blocks);

}

inline __host__ __device__ unsigned pow2(unsigned n) {
	return 1 << n;
}

}  // namespace asUtils


#endif /* HELPER_H_ */
