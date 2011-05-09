typedef double Rotmat[9];
typedef double Vec[3];

inline void matmult(Rotmat &rotmat1, Rotmat &rotmat2, Rotmat &rotmat) {
   Rotmat &m1 = rotmat1; Rotmat &m2 = rotmat2; Rotmat &m = rotmat;
   m[0] = m1[0]*m2[0]+m1[3]*m2[1]+m1[6]*m2[2];
   m[1] = m1[1]*m2[0]+m1[4]*m2[1]+m1[7]*m2[2];
   m[2] = m1[2]*m2[0]+m1[5]*m2[1]+m1[8]*m2[2];
   m[3] = m1[0]*m2[3]+m1[3]*m2[4]+m1[6]*m2[5];
   m[4] = m1[1]*m2[3]+m1[4]*m2[4]+m1[7]*m2[5];
   m[5] = m1[2]*m2[3]+m1[5]*m2[4]+m1[8]*m2[5];
   m[6] = m1[0]*m2[6]+m1[3]*m2[7]+m1[6]*m2[8];
   m[7] = m1[1]*m2[6]+m1[4]*m2[7]+m1[7]*m2[8];
   m[8] = m1[2]*m2[6]+m1[5]*m2[7]+m1[8]*m2[8];
}

extern "C" void matmult_(double *rotmat1, double *rotmat2, double *rotmat) {
  matmult(
    *((Rotmat *) rotmat1),
    *((Rotmat *) rotmat2),
    *((Rotmat *) rotmat)
  );
}

inline void vecmatmult(Vec &v0, Rotmat &m, Vec &v) {
  v[0] = v0[0]*m[0]+v0[1]*m[1]+v0[2]*m[2];
  v[1] = v0[0]*m[3]+v0[1]*m[4]+v0[2]*m[5];
  v[2] = v0[0]*m[6]+v0[1]*m[7]+v0[2]*m[8];
}

extern "C" void vecmatmult_(double *v0, double *m, double *v) {
  vecmatmult(
    *((Vec *) v0),
    *((Rotmat *) m),
    *((Vec *) v)
  );
}
