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

#ifndef UTILS_H_
#define UTILS_H_

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cuda_runtime.h>


#ifndef NDEBUG
#define cudaVerify(x) do { 																				\
		cudaError_t __cu_result = x; 																	\
		if (__cu_result!=cudaSuccess) { 																\
			fprintf(stderr, "%s:%i: Error: cuda function call failed:\n" 								\
					"%s;\nmessage: %s\n", 																\
					__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));							\
			exit(1);																					\
		} 																								\
	} while(0)
#define cudaVerifyKernel(x) do {																		\
		x;																								\
		cudaError_t __cu_result = cudaGetLastError();													\
		if (__cu_result!=cudaSuccess) { 																\
			fprintf(stderr, "%s:%i: Error: cuda function call failed:\n" 								\
					"%s;\nmessage: %s\n", 																\
					__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));							\
			exit(1);																					\
		} 																								\
	} while(0)
#else
#define cudaVerify(x) do {																				\
		x;																								\
	} while(0)
#define cudaVerifyKernel(x) do {																		\
		x;																								\
	} while(0)
#endif

#define CUDA_CHECK(x) do {                                                  \
	x;                                                                      \
	cudaError_t __cu_result = cudaGetLastError();							\
	if (__cu_result!=cudaSuccess) { 										\
		fprintf(stderr, "%s:%i: Error: cuda function call failed:\n" 		\
				"%s;\nmessage: %s\n", 										\
				__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));	\
		exit(1);															\
	} 																		\
} while(0)

#define ASSERT(x) do {                                           \
	if((x) == false) {                                           \
		fprintf(stderr, "%s:%i: Assertion '%s' failed.\n",     	\
				__FILE__, __LINE__, #x );                        \
		exit(1);                                                 \
	}                                                            \
} while(0)


#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) ((x < y) ? y : x)
#endif


#endif /* UTILS_H_ */
