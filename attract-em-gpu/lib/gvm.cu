#include "stdio.h"

typedef float Map[{{dimensions[0]}}][{{dimensions[1]}}][{{dimensions[2]}}];
typedef float RefeMap[{{dimensions2[0]}}][{{dimensions2[1]}}][{{dimensions2[2]}}];

__global__ void gvm(Map *maps, RefeMap *refe_map0, float *sumx, float *sumxx, float *sumxy) {
  if (threadIdx.x == 0 && blockIdx.y == 0) {    
    sumx[blockIdx.x] = 0;    
  }
  else {
    if (threadIdx.x == 1 && blockIdx.y == 0) {    
      sumxx[blockIdx.x] = 0;    
    }
    else if (threadIdx.x == 2 && blockIdx.y == 0) {    
      sumxy[blockIdx.x] = 0;    
    }
  }
  RefeMap &refe_map = *refe_map0;
  int blockz = blockIdx.y / {{dimensions[0]-2}};
  int blockx = blockIdx.y % {{dimensions[0]-2}};
  Map &map = maps[blockIdx.x];
  int offset_z = 30 * blockz;
  int offset_x = blockx;
  
  float plus[3][3];
  float minus[3][3];
  float delta = 0, refe = 0, xy = 0;
  
  int z = offset_z + threadIdx.x;
  
  if (z >= {{dimensions[blockdim]-2}}) return;
    
  for (int d1 = 0; d1 < 3; d1++) {
    int x = offset_x;
    int y = 0;
    for (int d2 = 0; d2 < 3; d2++) {
      minus[d1][d2] = map[{{DX}}][{{DY}}][{{DZ}}];
      plus[d1][d2] = map[{{DX}}{{TX}}][{{DY}}{{TY}}][{{DZ}}{{TZ}}]; 
    }
  }  

  for (int offset_y = 0; offset_y < {{dimensions[propdim] - 2}}; offset_y++) {
    refe = refe_map{{refemapcode}};
        
    float dplus = 0;
    float dminus = 0;
    for (int d1 = 0; d1 < 3; d1++) {
      for (int d2 = 0; d2 < 3; d2++) {
        dplus += plus[d1][d2];
        dminus += minus[d1][d2];
      }
    }

    // Calculate correlation part
    delta = 0;
    xy = 0;
    if (refe != 0) { //do this AFTER dplus and dminus have been computed, fetching refe is slow...
      delta = dplus - dminus;
      xy = delta * refe;
    }
    float sum_delta = delta, sum_xx = delta*delta, sum_xy = xy;
    for (int offset = 16; offset > 0; offset /= 2) {
      sum_delta += __shfl_down(sum_delta, offset);
      sum_xx += __shfl_down(sum_xx, offset);
      sum_xy += __shfl_down(sum_xy, offset);
    }    
    if (threadIdx.x == 0) {    
      atomicAdd(&sumx[blockIdx.x], sum_delta);
      atomicAdd(&sumxx[blockIdx.x], sum_xx);
      atomicAdd(&sumxy[blockIdx.x], sum_xy);      
    }

    // Fetch next data
    if (offset_y == {{dimensions[1] - 3}}) break;
    
    for (int d1 = 0; d1 < 3; d1++) {
      for (int d2 = 0; d2 < 2; d2++) {
        minus[d1][d2] = minus[d1][d2+1];
        plus[d1][d2] = plus[d1][d2+1];
      }
    }  
    for (int d1 = 0; d1 < 3; d1++) {      
      int x = offset_x;
      int y = offset_y + 1;
      int d2 = 2;      
      minus[d1][2] = map[{{DX}}][{{DY}}][{{DZ}}];
      plus[d1][2] = map[{{DX}}{{TX}}][{{DY}}{{TY}}][{{DZ}}{{TZ}}];
    }

    
  }
}  
