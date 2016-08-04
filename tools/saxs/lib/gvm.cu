#include "stdio.h"

typedef float Map[{{dimensions[0]}}][{{dimensions[1]}}][{{dimensions[2]}}];
typedef float RefeMap[{{dimensions2[0]}}][{{dimensions2[1]}}][{{dimensions2[2]}}];

__global__ void gvm_x(Map *maps, RefeMap *refe_map0, float *sumx, float *sumxx, float *sumxy) {
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
  int blocky = blockIdx.y / {{dimensions[2]-2}};
  int blockz = blockIdx.y % {{dimensions[2]-2}};
  Map &map = maps[blockIdx.x];
  int offset_x = 30 * blocky;
  int offset_z = blockz;
  
  float plus[3][3];
  float minus[3][3];
  float delta = 0, refe = 0, xy = 0;
  
  int x = offset_x + threadIdx.x;
  
  if (x + 2 >= {{dimensions[0]}}) return;
  
  for (int dy = 0; dy < 3; dy++) {
    for (int dz = 0; dz < 3; dz++) {
      int z = offset_z + dz;
      minus[dy][dz] = map[x][dy][z];
      plus[dy][dz] = map[x+2][dy][z]; 
    }
  }  

  for (int offset_y = 0; offset_y < {{dimensions[1] - 2}}; offset_y++) {
    refe = refe_map[x][offset_y][offset_z];
        
    float dplus = 0;
    float dminus = 0;
    for (int dy = 0; dy < 3; dy++) {
      for (int dz = 0; dz < 3; dz++) {
        dplus += plus[dy][dz];
        dminus += minus[dy][dz];
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
    
    for (int dy = 0; dy < 2; dy++) {
      for (int dz = 0; dz < 3; dz++) {
        minus[dy][dz] = minus[dy+1][dz];
        plus[dy][dz] = plus[dy+1][dz];
      }
    }  
    int y = offset_y + 3;
    for (int dz = 0; dz < 3; dz++) {
      int z = offset_z + dz;
      minus[2][dz] = map[x][y][z];
      plus[2][dz] = map[x+2][y][z];
    }

    
  }
}  

__global__ void gvm_y(Map *maps, RefeMap *refe_map0, float *sumx, float *sumxx, float *sumxy) {
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
  int blocky = blockIdx.y / {{dimensions[1]-2}};
  int blockz = blockIdx.y % {{dimensions[1]-2}};
  Map &map = maps[blockIdx.x];
  int offset_x = 30 * blocky + threadIdx.x;
  int y = blockz;
  
  float plus[3][3];
  float minus[3][3];
  float delta = 0, refe = 0, xy = 0;
  
  if (offset_x + 2 >= {{dimensions[0]}}) return;
  
  for (int dx = 0; dx < 3; dx++) {
    int x = offset_x + dx;
    for (int dz = 0; dz < 3; dz++) {      
      minus[dx][dz] = map[x][y][dz];
      plus[dx][dz] = map[x][y+2][dz]; 
    }
  }  

  for (int offset_z = 0; offset_z < {{dimensions[2] - 2}}; offset_z++) {
    refe = refe_map[offset_x][y][offset_z];
        
    float dplus = 0;
    float dminus = 0;
    for (int dx = 0; dx < 3; dx++) {
      for (int dz = 0; dz < 3; dz++) {
        dplus += plus[dx][dz];
        dminus += minus[dx][dz];
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
    if (offset_z == {{dimensions[2] - 3}}) break;
    
    for (int dz = 0; dz < 2; dz++) {
      for (int dx = 0; dx < 3; dx++) {
        minus[dx][dz] = minus[dx][dz+1];
        plus[dx][dz] = plus[dx][dz+1];
      }
    }  
    int z = offset_z + 3;
    for (int dx = 0; dx < 3; dx++) {
      int x = offset_x + dx;
      minus[dx][2] = map[x][y][z];
      plus[dx][2] = map[x][y+2][z];
    }

    
  }
}  

__global__ void gvm_z(Map *maps, RefeMap *refe_map0, float *sumx, float *sumxx, float *sumxy) {
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
  int blocky = blockIdx.y / {{dimensions[1]-2}};
  int blockz = blockIdx.y % {{dimensions[1]-2}};
  Map &map = maps[blockIdx.x];
  int offset_x = 30 * blocky + threadIdx.x;
  int z = blockz;
  
  float plus[3][3];
  float minus[3][3];
  float delta = 0, refe = 0, xy = 0;
  
  if (offset_x + 2 >= {{dimensions[0]}}) return;
  
  for (int dx = 0; dx < 3; dx++) {
    int x = offset_x + dx;
    for (int dy = 0; dy < 3; dy++) {      
      minus[dx][dy] = map[x][dy][z];
      plus[dx][dy] = map[x][dy][z+2]; 
    }
  }  

  for (int offset_y = 0; offset_y < {{dimensions[1] - 2}}; offset_y++) {
    refe = refe_map[offset_x][offset_y][z];
    
    // Calculate dplus and dminus
    float dplus = 0;
    float dminus = 0;
    for (int dx = 0; dx < 3; dx++) {
      for (int dy = 0; dy < 3; dy++) {
        dplus += plus[dx][dy];
        dminus += minus[dx][dy];
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
    
    for (int dy = 0; dy < 2; dy++) {
      for (int dx = 0; dx < 3; dx++) {
        minus[dx][dy] = minus[dx][dy+1];
        plus[dx][dy] = plus[dx][dy+1];
      }
    }  
    int y = offset_y + 3;
    for (int dx = 0; dx < 3; dx++) {
      int x = offset_x + dx;
      minus[dx][2] = map[x][y][z];
      plus[dx][2] = map[x][y][z+2];
    }

    
  }
}  
