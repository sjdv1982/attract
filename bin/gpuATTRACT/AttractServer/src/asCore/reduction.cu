
#include <stdio.h>

#include "asCore/reduction.h"
#include "as/asTypes.h"
#include "asUtils/helper.h"
#include "config.h"

extern __constant__ as::deviceProteinDesc c_Proteins[DEVICE_MAXPROTEINS];

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};


template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce6(T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl, T* d_x, T* d_y, T* d_z,  T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x; // half the number of blocks
    unsigned int gridSize = blockSize*2*gridDim.x;

    T sum_fx = 0;
    T sum_fy = 0;
    T sum_fz = 0;
    T sum_eVdW = 0;
    T sum_eEl = 0;
    T sum_torque[9] = {0};

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
    	T fx,fy,fz, x, y, z;
    	fx = d_fx[i];
    	fy = d_fy[i];
    	fz = d_fz[i];
    	x = d_x[i];
    	y = d_y[i];
    	z = d_z[i];
        sum_fx 	 += fx;
        sum_fy 	 += fy;
        sum_fz 	 += fz;
        sum_eVdW 	 += d_eVdW[i];
        sum_eEl 	 += d_eEl[i];
        sum_torque[0] += x*fx;
		sum_torque[1] += y*fx;
		sum_torque[2] += z*fx;
		sum_torque[3] += x*fy;
		sum_torque[4] += y*fy;
		sum_torque[5] += z*fy;
		sum_torque[6] += x*fz;
		sum_torque[7] += y*fz;
		sum_torque[8] += z*fz;


        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n) {

			fx = d_fx[i+blockSize];
			fy = d_fy[i+blockSize];
			fz = d_fz[i+blockSize];
			x = d_x[i+blockSize];
			y = d_y[i+blockSize];
			z = d_z[i+blockSize];
			sum_fx 	 += fx;
			sum_fy 	 += fy;
			sum_fz 	 += fz;
			sum_eVdW 	 += d_eVdW[i+blockSize];
			sum_eEl 	 += d_eEl[i+blockSize];
			sum_torque[0] += x*fx;
			sum_torque[1] += y*fx;
			sum_torque[2] += z*fx;
			sum_torque[3] += x*fy;
			sum_torque[4] += y*fy;
			sum_torque[5] += z*fy;
			sum_torque[6] += x*fz;
			sum_torque[7] += y*fz;
			sum_torque[8] += z*fz;
        }


        i += gridSize;
    }

    // each thread puts its local sum into shared memory
    sdata[tid + 0 * blockSize] = sum_fx;
    sdata[tid + 1 * blockSize] = sum_fy;
    sdata[tid + 2 * blockSize] = sum_fz;
    sdata[tid + 3 * blockSize] = sum_eVdW;
    sdata[tid + 4 * blockSize] = sum_eEl;
    sdata[tid + 5 * blockSize] = sum_torque[0];
    sdata[tid + 6 * blockSize] = sum_torque[1];
    sdata[tid + 7 * blockSize] = sum_torque[2];
    sdata[tid + 8 * blockSize] = sum_torque[3];
    sdata[tid + 9 * blockSize] = sum_torque[4];
    sdata[tid + 10* blockSize] = sum_torque[5];
    sdata[tid + 11* blockSize] = sum_torque[6];
    sdata[tid + 12* blockSize] = sum_torque[7];
    sdata[tid + 13* blockSize] = sum_torque[8];

    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 1024) && (tid < 512))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 512];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 512];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 512];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 512];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 512];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 512];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 512];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 512];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 512];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 512];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 512];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 512];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 512];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 512];
    }
    
    __syncthreads();

    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 256];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 256];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 256];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 256];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 256];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 256];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 256];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 256];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 256];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 256];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 256];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 256];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 256];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 256];
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
		sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 128];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 128];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 128];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 128];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 128];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 128];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 128];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 128];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 128];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 128];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 128];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 128];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 128];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 128];
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 64];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 64];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 64];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 64];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 64];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 64];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 64];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 64];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 64];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 64];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 64];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 64];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 64];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 64];
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) {
        	sum_fx        += sdata[tid + 0 * blockSize + 32];
			sum_fy        += sdata[tid + 1 * blockSize + 32];
			sum_fz        += sdata[tid + 2 * blockSize + 32];
			sum_eVdW      += sdata[tid + 3 * blockSize + 32];
			sum_eEl       += sdata[tid + 4 * blockSize + 32];
			sum_torque[0] += sdata[tid + 5 * blockSize + 32];
			sum_torque[1] += sdata[tid + 6 * blockSize + 32];
			sum_torque[2] += sdata[tid + 7 * blockSize + 32];
			sum_torque[3] += sdata[tid + 8 * blockSize + 32];
			sum_torque[4] += sdata[tid + 9 * blockSize + 32];
			sum_torque[5] += sdata[tid + 10* blockSize + 32];
			sum_torque[6] += sdata[tid + 11* blockSize + 32];
			sum_torque[7] += sdata[tid + 12* blockSize + 32];
			sum_torque[8] += sdata[tid + 13* blockSize + 32];

        }
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            sum_fx        += __shfl_down(sum_fx       , offset);
			sum_fy        += __shfl_down(sum_fy       , offset);
			sum_fz        += __shfl_down(sum_fz       , offset);
			sum_eVdW      += __shfl_down(sum_eVdW     , offset);
			sum_eEl       += __shfl_down(sum_eEl      , offset);
			sum_torque[0] += __shfl_down(sum_torque[0], offset);
			sum_torque[1] += __shfl_down(sum_torque[1], offset);
			sum_torque[2] += __shfl_down(sum_torque[2], offset);
			sum_torque[3] += __shfl_down(sum_torque[3], offset);
			sum_torque[4] += __shfl_down(sum_torque[4], offset);
			sum_torque[5] += __shfl_down(sum_torque[5], offset);
			sum_torque[6] += __shfl_down(sum_torque[6], offset);
			sum_torque[7] += __shfl_down(sum_torque[7], offset);
			sum_torque[8] += __shfl_down(sum_torque[8], offset);
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 32];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 32];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 32];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 32];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 32];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 32];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 32];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 32];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 32];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 32];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 32];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 32];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 32];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 32];
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 16];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 16];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 16];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 16];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 16];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 16];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 16];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 16];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 16];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 16];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 16];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 16];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 16];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 16];
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 8];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 8];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 8];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 8];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 8];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 8];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 8];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 8];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 8];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 8];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 8];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 8];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 8];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 8];
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 4];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 4];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 4];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 4];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 4];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 4];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 4];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 4];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 4];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 4];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 4];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 4];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 4];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 4];
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 2];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 2];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 2];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 2];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 2];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 2];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 2];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 2];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 2];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 2];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 2];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 2];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 2];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 2];
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 1];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 1];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 1];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 1];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 1];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 1];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 1];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 1];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 1];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 1];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 1];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 1];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 1];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 1];
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) {
    	g_odata[gridDim.x * 0  + blockIdx.x] = sum_fx;
    	g_odata[gridDim.x * 1  + blockIdx.x] = sum_fy;
    	g_odata[gridDim.x * 2  + blockIdx.x] = sum_fz;
    	g_odata[gridDim.x * 3  + blockIdx.x] = sum_eVdW;
    	g_odata[gridDim.x * 4  + blockIdx.x] = sum_eEl;
    	g_odata[gridDim.x * 5  + blockIdx.x ] = sum_torque[0];
    	g_odata[gridDim.x * 6  + blockIdx.x ] = sum_torque[1];
    	g_odata[gridDim.x * 7  + blockIdx.x ] = sum_torque[2];
    	g_odata[gridDim.x * 8  + blockIdx.x ] = sum_torque[3];
    	g_odata[gridDim.x * 9  + blockIdx.x ] = sum_torque[4];
    	g_odata[gridDim.x * 10 + blockIdx.x ] = sum_torque[5];
    	g_odata[gridDim.x * 11 + blockIdx.x ] = sum_torque[6];
    	g_odata[gridDim.x * 12 + blockIdx.x ] = sum_torque[7];
    	g_odata[gridDim.x * 13 + blockIdx.x ] = sum_torque[8];
    }
}


template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce1Grid(unsigned int protId,
		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl, T *g_odata)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;

    unsigned int i = threadIdx.x; // half the number of blocks
    unsigned int base = blockIdx.x*c_Proteins[protId].numAtoms;

    T sum_fx = 0;
    T sum_fy = 0;
    T sum_fz = 0;
    T sum_eVdW = 0;
    T sum_eEl = 0;
    T sum_torque[9] = {0};

    while (i < c_Proteins[protId].numAtoms)
    {
		T fx, fy, fz, x, y, z;
		fx = d_fx[base + i];
		fy = d_fy[base + i];
		fz = d_fz[base + i];
		x = c_Proteins[protId].xPos[i];
		y = c_Proteins[protId].yPos[i];
		z = c_Proteins[protId].zPos[i];
		sum_fx += fx;
		sum_fy += fy;
		sum_fz += fz;
		sum_eVdW += d_eVdW[base + i];
		sum_eEl += d_eEl[base + i];
		sum_torque[0] += x * fx;
		sum_torque[1] += y * fx;
		sum_torque[2] += z * fx;
		sum_torque[3] += x * fy;
		sum_torque[4] += y * fy;
		sum_torque[5] += z * fy;
		sum_torque[6] += x * fz;
		sum_torque[7] += y * fz;
		sum_torque[8] += z * fz;


        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < c_Proteins[protId].numAtoms) {

			fx = d_fx[base + i + blockSize];
			fy = d_fy[base + i + blockSize];
			fz = d_fz[base + i + blockSize];
			x = c_Proteins[protId].xPos[i + blockSize];
			y = c_Proteins[protId].yPos[i + blockSize];
			z = c_Proteins[protId].zPos[i + blockSize];
			sum_fx += fx;
			sum_fy += fy;
			sum_fz += fz;
			sum_eVdW += d_eVdW[base + i + blockSize];
			sum_eEl += d_eEl[base + i + blockSize];
			sum_torque[0] += x * fx;
			sum_torque[1] += y * fx;
			sum_torque[2] += z * fx;
			sum_torque[3] += x * fy;
			sum_torque[4] += y * fy;
			sum_torque[5] += z * fy;
			sum_torque[6] += x * fz;
			sum_torque[7] += y * fz;
			sum_torque[8] += z * fz;
		}


        i += blockSize*2;
    }



    // each thread puts its local sum into shared memory
    sdata[tid + 0 * blockSize] = sum_fx;
    sdata[tid + 1 * blockSize] = sum_fy;
    sdata[tid + 2 * blockSize] = sum_fz;
    sdata[tid + 3 * blockSize] = sum_eVdW;
    sdata[tid + 4 * blockSize] = sum_eEl;
    sdata[tid + 5 * blockSize] = sum_torque[0];
    sdata[tid + 6 * blockSize] = sum_torque[1];
    sdata[tid + 7 * blockSize] = sum_torque[2];
    sdata[tid + 8 * blockSize] = sum_torque[3];
    sdata[tid + 9 * blockSize] = sum_torque[4];
    sdata[tid + 10* blockSize] = sum_torque[5];
    sdata[tid + 11* blockSize] = sum_torque[6];
    sdata[tid + 12* blockSize] = sum_torque[7];
    sdata[tid + 13* blockSize] = sum_torque[8];

    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 1024) && (tid < 512))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 512];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 512];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 512];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 512];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 512];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 512];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 512];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 512];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 512];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 512];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 512];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 512];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 512];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 512];
    }
    
	__syncthreads();
	
    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 256];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 256];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 256];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 256];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 256];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 256];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 256];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 256];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 256];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 256];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 256];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 256];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 256];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 256];
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
		sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 128];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 128];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 128];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 128];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 128];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 128];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 128];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 128];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 128];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 128];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 128];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 128];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 128];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 128];
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 64];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 64];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 64];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 64];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 64];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 64];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 64];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 64];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 64];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 64];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 64];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 64];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 64];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 64];
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) {
        	sum_fx        += sdata[tid + 0 * blockSize + 32];
			sum_fy        += sdata[tid + 1 * blockSize + 32];
			sum_fz        += sdata[tid + 2 * blockSize + 32];
			sum_eVdW      += sdata[tid + 3 * blockSize + 32];
			sum_eEl       += sdata[tid + 4 * blockSize + 32];
			sum_torque[0] += sdata[tid + 5 * blockSize + 32];
			sum_torque[1] += sdata[tid + 6 * blockSize + 32];
			sum_torque[2] += sdata[tid + 7 * blockSize + 32];
			sum_torque[3] += sdata[tid + 8 * blockSize + 32];
			sum_torque[4] += sdata[tid + 9 * blockSize + 32];
			sum_torque[5] += sdata[tid + 10* blockSize + 32];
			sum_torque[6] += sdata[tid + 11* blockSize + 32];
			sum_torque[7] += sdata[tid + 12* blockSize + 32];
			sum_torque[8] += sdata[tid + 13* blockSize + 32];

        }
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            sum_fx        += __shfl_down(sum_fx       , offset);
			sum_fy        += __shfl_down(sum_fy       , offset);
			sum_fz        += __shfl_down(sum_fz       , offset);
			sum_eVdW      += __shfl_down(sum_eVdW     , offset);
			sum_eEl       += __shfl_down(sum_eEl      , offset);
			sum_torque[0] += __shfl_down(sum_torque[0], offset);
			sum_torque[1] += __shfl_down(sum_torque[1], offset);
			sum_torque[2] += __shfl_down(sum_torque[2], offset);
			sum_torque[3] += __shfl_down(sum_torque[3], offset);
			sum_torque[4] += __shfl_down(sum_torque[4], offset);
			sum_torque[5] += __shfl_down(sum_torque[5], offset);
			sum_torque[6] += __shfl_down(sum_torque[6], offset);
			sum_torque[7] += __shfl_down(sum_torque[7], offset);
			sum_torque[8] += __shfl_down(sum_torque[8], offset);
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 32];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 32];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 32];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 32];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 32];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 32];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 32];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 32];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 32];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 32];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 32];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 32];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 32];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 32];
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 16];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 16];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 16];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 16];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 16];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 16];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 16];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 16];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 16];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 16];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 16];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 16];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 16];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 16];
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 8];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 8];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 8];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 8];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 8];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 8];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 8];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 8];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 8];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 8];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 8];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 8];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 8];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 8];
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 4];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 4];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 4];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 4];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 4];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 4];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 4];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 4];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 4];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 4];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 4];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 4];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 4];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 4];
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 2];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 2];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 2];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 2];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 2];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 2];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 2];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 2];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 2];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 2];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 2];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 2];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 2];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 2];
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 1];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 1];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 1];
		sdata[tid + 3 * blockSize] = sum_eVdW      = sum_eVdW      + sdata[tid + 3 * blockSize + 1];
		sdata[tid + 4 * blockSize] = sum_eEl       = sum_eEl       + sdata[tid + 4 * blockSize + 1];
		sdata[tid + 5 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 5 * blockSize + 1];
		sdata[tid + 6 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 6 * blockSize + 1];
		sdata[tid + 7 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 7 * blockSize + 1];
		sdata[tid + 8 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 8 * blockSize + 1];
		sdata[tid + 9 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 9 * blockSize + 1];
		sdata[tid + 10* blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 10* blockSize + 1];
		sdata[tid + 11* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 11* blockSize + 1];
		sdata[tid + 12* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 12* blockSize + 1];
		sdata[tid + 13* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 13* blockSize + 1];
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) {
    	g_odata[0  + blockIdx.x*14] = sum_fx;
    	g_odata[1  + blockIdx.x*14] = sum_fy;
    	g_odata[2  + blockIdx.x*14] = sum_fz;
    	g_odata[3  + blockIdx.x*14] = sum_eVdW;
    	g_odata[4  + blockIdx.x*14] = sum_eEl;
    	g_odata[5  + blockIdx.x*14] = sum_torque[0];
    	g_odata[6  + blockIdx.x*14] = sum_torque[1];
    	g_odata[7  + blockIdx.x*14] = sum_torque[2];
    	g_odata[8  + blockIdx.x*14] = sum_torque[3];
    	g_odata[9  + blockIdx.x*14] = sum_torque[4];
    	g_odata[10 + blockIdx.x*14] = sum_torque[5];
    	g_odata[11 + blockIdx.x*14] = sum_torque[6];
    	g_odata[12 + blockIdx.x*14] = sum_torque[7];
    	g_odata[13 + blockIdx.x*14] = sum_torque[8];
    }
}


////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////


//template <class T>
//inline __device__ void childReduce(int size, int threads, int blocks,
//		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
//		T* d_x, T* d_y, T* d_z, T *d_odata,
//		const cudaStream_t stream = 0)
//{
//	dim3 dimBlock(threads, 1, 1);
//	dim3 dimGrid(blocks, 1, 1);
//
//	// when there is only one warp per block, we need to allocate two warps
//	// worth of shared memory so that we don't index shared memory out of bounds
//	int smemSize = (threads <= 32) ? 2 * 14 * threads * sizeof(T) : 14 * threads * sizeof(T);
//
//	// choose which of the optimized versions of reduction to launch
//
//	if (asUtils::isPow2(size))
//	{
//		switch (threads)
//		{
//			case 1024:
//				reduce6<T, 1024,true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 512:
//				reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 256:
//				reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 128:
//				reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 64:
//				reduce6<T,  64, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 32:
//				reduce6<T,  32, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 16:
//				reduce6<T,  16, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  8:
//				reduce6<T,   8, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  4:
//				reduce6<T,   4, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  2:
//				reduce6<T,   2, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  1:
//				reduce6<T,   1, true><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//		}
//	}
//	else
//	{
//		switch (threads)
//		{
//			case 1024:
//				reduce6<T, 1024,false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 512:
//				reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 256:
//				reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 128:
//				reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 64:
//				reduce6<T,  64, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 32:
//				reduce6<T,  32, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case 16:
//				reduce6<T,  16, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  8:
//				reduce6<T,   8, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  4:
//				reduce6<T,   4, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  2:
//				reduce6<T,   2, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//
//			case  1:
//				reduce6<T,   1, false><<< dimGrid, dimBlock, smemSize, stream >>>(d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_x, d_y, d_z, d_odata, size);
//				break;
//		}
//	}
//
//
//}


//reduce1Grid(unsigned int protId, unsigned int n,
//		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl, T *g_odata)

template <class T>
__host__ void asCore::reduceAll(const unsigned& threads, const unsigned& blocks,
		const unsigned& protId, const unsigned& size,
		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
		T *d_odata,
		const cudaStream_t& stream)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * 14 * threads * sizeof(T) : 14 * threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch

	if (asUtils::isPow2(size))
	{
		switch (threads)
		{
			case 1024:
				reduce1Grid<T, 1024,true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 512:
				reduce1Grid<T, 512, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 256:
				reduce1Grid<T, 256, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 128:
				reduce1Grid<T, 128, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 64:
				reduce1Grid<T,  64, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 32:
				reduce1Grid<T,  32, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 16:
				reduce1Grid<T,  16, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  8:
				reduce1Grid<T,   8, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  4:
				reduce1Grid<T,   4, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  2:
				reduce1Grid<T,   2, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  1:
				reduce1Grid<T,   1, true><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;
		}
	}
	else
	{
		switch (threads)
		{
			case 1024:
				reduce1Grid<T, 1024,false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 512:
				reduce1Grid<T, 512, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 256:
				reduce1Grid<T, 256, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 128:
				reduce1Grid<T, 128, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 64:
				reduce1Grid<T,  64, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 32:
				reduce1Grid<T,  32, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case 16:
				reduce1Grid<T,  16, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  8:
				reduce1Grid<T,   8, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  4:
				reduce1Grid<T,   4, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  2:
				reduce1Grid<T,   2, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;

			case  1:
				reduce1Grid<T,   1, false><<< dimGrid, dimBlock, smemSize, stream >>>(protId, d_fx, d_fy, d_fz, d_eVdW, d_eEl, d_odata);
				break;
		}
	}


}

/*
 ** @remark: in order to support normal modes it is not enough to supply the protId.
 ** With the protId you can gather only undeformed coordinates. However, for reduction, deformed
 ** coordinates are required.
 ** To integrate normal modes into the pipeline means that one needs to separate deformation
 ** and rotation+translation in terms of execution and memory storage. Otherwise deformation
 ** needs to be calculated twice.
 */

//template <class T>
//__global__ void asCore::reduce(unsigned protId,
//       T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
//       T *d_odata)
//{
//	/* size: number of elements per object (numAtoms) */
//
//	/* The first thread in a block starts a new reduction for a specific structure.
//	 * The GridDim == chunkSize */
//
//	unsigned int tid = threadIdx.x;
//	if (tid == 0) {
//		const int size = c_Proteins[protId].numAtoms;
//
//		// 512 == max blockSize -> extension to 1024 would lead to too much shared memory requirements!!!!
//		// 1 max number of blocks. Do not change it!
//		int threads = 0;
//		int blocks = 0;
//		asUtils::getNumBlocksAndThreads(size, 1 , 512, blocks, threads);
//
//		childReduce<T>(size,threads, blocks,
//					d_fx 	+ blockIdx.x*size,
//					d_fy 	+ blockIdx.x*size,
//					d_fz 	+ blockIdx.x*size,
//					d_eVdW 	+ blockIdx.x*size,
//					d_eEl 	+ blockIdx.x*size,
//					c_Proteins[protId].xPos,
//					c_Proteins[protId].yPos,
//					c_Proteins[protId].zPos,
//					d_odata + blockIdx.x*14);
//
//#ifndef NDEBUG
//		if(blockIdx.x == 0) {
//			cudaDeviceSynchronize();
//			cudaError_t __cu_result = cudaGetLastError();
//			if (__cu_result!=cudaSuccess) {
//				printf("%s:%i: error: cuda function call failed:\n;\nmessage: %s\n",
//						__FILE__, __LINE__, cudaGetErrorString(__cu_result) );
//				return;
//			}
//		}
//#endif
//	}
//}




/* explicit instantiation */
//template __global__ void
//asCore::reduce<float>(unsigned protId,
//       float *d_fx, float *d_fy, float *d_fz, float *d_eVdW, float *d_eEl,
//       float *d_odata);

template __host__ void
asCore::reduceAll<float>(const unsigned& threads, const unsigned& blocks,
		const unsigned& protId, const unsigned&,
		float *d_fx, float *d_fy, float *d_fz, float *d_eVdW, float *d_eEl,
		float *d_odata,
		const cudaStream_t& stream);




