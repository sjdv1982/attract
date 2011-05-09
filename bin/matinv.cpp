#include <cstdlib>
#include <cstring>
using namespace std;

typedef double Rotmat[9];

inline void matinv(Rotmat &rotmat) {
  Rotmat m0;
  memcpy(m0, rotmat, sizeof(Rotmat));
  Rotmat &m = rotmat;
  //m[0] = m0[0];
  m[1] = m0[3];
  m[2] = m0[6];
  m[3] = m0[1];
  //m[4] = m0[4];
  m[5] = m0[7];
  m[6] = m0[2];
  m[7] = m0[5];
  //m[8] = m0[8];
}

extern "C" void matinv_(double *rotmat) {
  matinv (
    *((Rotmat *) rotmat)
  );
}

inline void matcopy(Rotmat &rotmat, Rotmat &rotmatcopy) {
  memcpy(rotmatcopy, rotmat, sizeof(Rotmat));
}

extern "C" void matcopy_(double *rotmat, double *rotmatcopy) {
  matcopy (
    *((Rotmat *) rotmat),
    *((Rotmat *) rotmatcopy)
  );
}

