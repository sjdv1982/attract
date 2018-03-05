#include "stdio.h"

typedef float Coordinate[3];

__device__ void trilin(float *grid, float weight, float ax, float ay, float az, int dim, float dimxx) {
  //dimxx = dim * dim
  float wx1=0,wy1=0,wz1=0;
    
  int px0 = floor(ax);  
  int px1 = ceil(ax);
  if (px1 < 0 || px0 >= dim) return;
  if (px0 < 0) wx1 = 1;
  else if (px1 >= dim) wx1 = 0;
  else wx1 = ax - px0;

  int py0 = floor(ay);  
  int py1 = ceil(ay);
  if (py1 < 0 || py0 >= dim) return;
  if (py0 < 0) wy1 = 1;
  else if (py1 >= dim) wy1 = 0;
  else wy1 = ay - py0;

  int pz0 = floor(az);  
  int pz1 = ceil(az);
  if (pz1 < 0 || pz0 >= dim) return;
  if (pz0 < 0) wz1 = 1;
  else if (pz1 >= dim) wz1 = 0;
  else wz1 = az - pz0;
  
  float wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;  
  float cwx, cwxy, cwxyz;
  int indx,indxy,indxyz; 
  if (wx0 > 0) {
    cwx = weight * wx0;
    indx = dimxx*px0;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dim*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;                   
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dim*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
  }
  if (wx1 > 0) {
    cwx = weight * wx1;
    indx = dimxx*px1;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dim*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);        
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dim*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        {{ debug }}
        atomicAdd(&grid[indxyz], cwxyz);
      }
    }
  }
}

__device__ float trilin2(float *grid, float ax, float ay, float az, int dim, float dimxx) {
  //dimxx = dim * dim
  float o = 0;
  float wx1=0,wy1=0,wz1=0;
    
  int px0 = floor(ax);  
  int px1 = ceil(ax);
  if (px1 < 0 || px0 >= dim) return 0;
  if (px0 < 0) wx1 = 1;
  else if (px1 >= dim) wx1 = 0;
  else wx1 = ax - px0;

  int py0 = floor(ay);  
  int py1 = ceil(ay);
  if (py1 < 0 || py0 >= dim) return 0;
  if (py0 < 0) wy1 = 1;
  else if (py1 >= dim) wy1 = 0;
  else wy1 = ay - py0;

  int pz0 = floor(az);  
  int pz1 = ceil(az);
  if (pz1 < 0 || pz0 >= dim) return 0;
  if (pz0 < 0) wz1 = 1;
  else if (pz1 >= dim) wz1 = 0;
  else wz1 = az - pz0;
  
  float wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;  
  float cwx, cwxy, cwxyz;
  int indx,indxy,indxyz; 
  float v;
  if (wx0 > 0) {
    cwx = wx0;
    indx = dimxx*px0;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dim*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;                   
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dim*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
    }
  }
  if (wx1 > 0) {
    cwx = wx1;
    indx = dimxx*px1;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + dim*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + dim*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + pz0;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + pz1;      
        v = grid[indxyz];
        if (v > 0) {
          {{ debug2 }}
          o += cwxyz*v*v/{{ maxdensity**2 }};          
        }  
      }
    }
  }
  return o;
}

__global__ void atomdensitymask1(
 Coordinate *coor, float *weights, int size_coor, 
 float *grids
) {  
  const int struc = blockIdx.x;
  const int c = blockIdx.y * blockDim.x + threadIdx.x;
  if (c >= size_coor) return;   
  const int grid_offset = struc * {{ dimension**3 }};
  const int offset = struc * size_coor + c;
  float ax = (coor[offset][0] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  float ay = (coor[offset][1] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  float az = (coor[offset][2] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  float weight = weights[c];
  trilin(grids+grid_offset, weight, ax, ay, az, {{ dimension }}, {{ dimension**2 }});
}

__global__ void atomdensitymask2(
 Coordinate *coor, float *overlaps, int size_coor, 
 float *grids
) {  
  const int struc = blockIdx.x;
  const int c = blockIdx.y * blockDim.x + threadIdx.x;

  if (c >= size_coor) return;   
  const int grid_offset = struc * {{ dimension**3 }};
  const int offset = struc * size_coor + c;
  float ax = (coor[offset][0] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  float ay = (coor[offset][1] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  float az = (coor[offset][2] + {{ 0.5 * dimension * gridspacing }}) / {{ gridspacing }};
  

  float o = trilin2(grids+grid_offset, ax, ay, az, {{ dimension }}, {{ dimension**2 }});
  overlaps[offset] = o;
}


