#include "stdio.h"

typedef float Coordinate[3];

const unsigned int coorlen1 = {{coorlen[0]}};
const unsigned int coorlen2 = {{coorlen[1]}};
const unsigned int nstruc = {{nstruc}};
const unsigned int nbins = {{nbins}};
const float binsize = {{binsize}};

const int THREADS_PER_BLOCK = {{THREADS_PER_BLOCK}};

typedef Coordinate Coordinate1[coorlen1];
typedef Coordinate1 Coordinates1[nstruc];
typedef Coordinate Coordinate2[coorlen2];
typedef Coordinate2 Coordinates2[nstruc];

typedef float Radist[nbins];
typedef Radist Radists[nstruc];

inline __device__ void calc_single_radist(
  Coordinate &c1, float w1,
  Coordinate &c2, float w2, 
  Radist &radist
) {
  Coordinate d;
  Coordinate dd;
  d[0] = (c1[0]-c2[0]);
  d[1] = (c1[1]-c2[1]);
  d[2] = (c1[2]-c2[2]);
  dd[0] = d[0]*d[0];
  dd[1] = d[1]*d[1];
  dd[2] = d[2]*d[2];
  float dissq = dd[0] + dd[1] + dd[2];
  //printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f  %.3f %.3f\n", dissq, c1[0], c1[1], c1[2], c2[0], c2[1], c2[2], w1, w2);
  unsigned int bin = dissq / binsize;
  if (bin < nbins) {
    //printf("%d %.3f\n", bin, w1*w2);
    atomicAdd(&radist[bin],  w1 * w2);
  }  
}

__global__ void calc_radist(
 Coordinates1 *coors1, float *weights1,
 Coordinates2 *coors2, float *weights2,
 Radists *radist
) {  
  const unsigned int struc = blockIdx.x;
  Coordinate1 &c1 = (*coors1)[struc];
  Coordinate2 &c2 = (*coors2)[struc];
  Radist &rad = (*radist)[struc];
  
  const unsigned int x_offset = THREADS_PER_BLOCK * blockIdx.y + threadIdx.x;
  float wr;
  Coordinate r;
  if (x_offset < coorlen1) {
    Coordinate &r0 = c1[x_offset];    
    r[0] = r0[0];
    r[1] = r0[1];
    r[2] = r0[2];
    wr = weights1[x_offset];
  }
  
  __shared__ Coordinate  sh_l[THREADS_PER_BLOCK];
  __shared__ float  sh_wl[THREADS_PER_BLOCK];
  for (unsigned int y = 0; y < coorlen2; y+= THREADS_PER_BLOCK) {
    
    __syncthreads();  
    
    int y_offset = y + threadIdx.x;
    if (y_offset < coorlen2) {
      Coordinate &l0 = c2[y_offset];
      Coordinate &l = sh_l[threadIdx.x];
      l[0] = l0[0];
      l[1] = l0[1];
      l[2] = l0[2];
      sh_wl[threadIdx.x] = weights2[y_offset];
    }
  
    __syncthreads();  
    
    if (x_offset < coorlen1) {
      unsigned int warpOffset = threadIdx.x / warpSize;
      for (unsigned int py0 = 0; py0 < THREADS_PER_BLOCK; py0++) {
        int py = (py0 + warpOffset) % THREADS_PER_BLOCK;
        if (y + py < coorlen2) {
          Coordinate &l = sh_l[py]; 
          float wl = sh_wl[py];
          calc_single_radist(r, wr, l, wl, rad);
        }  
      }
    }  
    
  }  
  
  
  
  
}