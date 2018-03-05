typedef float Euler[3];
typedef float RotMat[9];

__device__ void _euler2rotmat(float phi,float ssi,float rot, RotMat r) {
  float cs=cos(ssi);
  float cp=cos(phi);
  float ss=sin(ssi);
  float sp=sin(phi);
  float cscp=cs*cp;
  float cssp=cs*sp;
  float sscp=ss*cp;
  float sssp=ss*sp;
  float crot=cos(rot);
  float srot=sin(rot);

  r[0] = crot * cscp + srot * sp;
  r[1] = srot * cscp - crot * sp;
  r[2] = sscp;

  r[3] = crot * cssp - srot * cp;
  r[4] = srot * cssp + crot * cp;
  r[5] = sssp;

  r[6] = -crot * ss;
  r[7] = -srot * ss;
  r[8] = cs;
}  

__global__ void euler2rotmat(Euler *eulers, RotMat *rotmats, int arraylength) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < arraylength) {
    float phi = eulers[i][0];
    float ssi = eulers[i][1];
    float rot = eulers[i][2];
    RotMat &r = rotmats[i];
    _euler2rotmat(phi, ssi, rot, r);
  }
}