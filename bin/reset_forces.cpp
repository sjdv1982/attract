#include <cstdlib>
#include <cstring>
using namespace std;

extern "C" void cartstate_get_forces_(const int &handle,double *&f, int &nall3);

extern "C" void reset_forces_(const int &cartstatehandle) {
  double *f; int nall3;
  cartstate_get_forces_(cartstatehandle, f, nall3);
  memset(f, 0, nall3*sizeof(double));      
}
