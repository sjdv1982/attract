#include "stdio.h"

typedef float Coordinate[3];

__device__ void trilin(float *grid, float ax, float ay, float az, int dimx, int dimy, int dimz, float dimxy) {
  //dimxy = dimx * dimy
  float wx1=0,wy1=0,wz1=0;
    
  int px0 = floor(ax);  
  int px1 = ceil(ax);
  if (px1 < 0 || px0 >= dimx) return;
  if (px0 < 0) wx1 = 1;
  else if (px1 >= dimx) wx1 = 0;
  else wx1 = ax - px0;

  int py0 = floor(ay);  
  int py1 = ceil(ay);
  if (py1 < 0 || py0 >= dimy) return;
  if (py0 < 0) wy1 = 1;
  else if (py1 >= dimy) wy1 = 0;
  else wy1 = ay - py0;

  int pz0 = floor(az);  
  int pz1 = ceil(az);
  if (pz1 < 0 || pz0 >= dimz) return;
  if (pz0 < 0) wz1 = 1;
  else if (pz1 >= dimz) wz1 = 0;
  else wz1 = az - pz0;
  
  float wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;  
  float cwx, cwxy, cwxyz;
  int indx,indxy,indxyz; 
  if (wx0 > 0) {
    cwx = wx0;
    indx = px0;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;                   
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
  }
  if (wx1 > 0) {
    cwx = wx1;
    indx = px1;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);        
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
  }
}

__global__ void gridify(
 Coordinate *coor, int size_coor, 
 float *grids
) {  
  const int struc = blockIdx.x;
  const int c = blockIdx.y * blockDim.x + threadIdx.x;
  if (c >= size_coor) return;   
  const int grid_offset = struc * {{ dim[0] * dim[1] * dim[2] }};
  const int offset = struc * size_coor + c;
  float ax = (coor[offset][0] + {{ -origin[0] }} ) / {{ gridspacing }};
  float ay = (coor[offset][1] + {{ -origin[1] }}) / {{ gridspacing }};
  float az = (coor[offset][2] + {{ -origin[2] }}) / {{ gridspacing }};
  trilin(grids+grid_offset, az, ay, ax, {{ dim[2] }}, {{ dim[1] }}, {{ dim[0] }}, {{ dim[1] * dim[2] }});
}