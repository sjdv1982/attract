#include <ctime>
#include <cstdlib>

extern "C" void crand_(const int &seed, const int &size, double *r) {
  srand(seed);
  for (int n = 0; n < size; n++) {
    r[n] = double(rand())/RAND_MAX;
  }
}
