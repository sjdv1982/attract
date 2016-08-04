typedef float RotMat[9];
typedef float Coordinate[3];

__global__ void collect(
 RotMat *rotmats, 
 Coordinate *translations,  
 Coordinate *tmplate, 
 int size_template, 
 Coordinate *coordinates
) { 
  
  //Load rotation and translation matrices into shared memory
  __shared__ float m[9];
  __shared__ float c[3];
  const int struc = blockIdx.x;
  if (threadIdx.x < 12) {
    if (threadIdx.x < 9) {
      m[threadIdx.x] = rotmats[struc][threadIdx.x];
    }
    else {
      c[threadIdx.x-9] = translations[struc][threadIdx.x-9];
    }
  }
  __syncthreads();    

  int coor = blockIdx.y * blockDim.x + threadIdx.x;
  if (coor >= size_template) return;
  
  
  int offset = struc * size_template + coor;
  float x = tmplate[coor][0];
  float y = tmplate[coor][1];
  float z = tmplate[coor][2];
   
  coordinates[offset][0] = x*m[0]+y*m[1]+z*m[2]+c[0];
  coordinates[offset][1] = x*m[3]+y*m[4]+z*m[5]+c[1];
  coordinates[offset][2] = x*m[6]+y*m[7]+z*m[8]+c[2];  
  
}