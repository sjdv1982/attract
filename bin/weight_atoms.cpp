#include "cstdio"
#include "cstdlib"
using namespace std;

const double aweights[] = {
 0, //0
 27.0, //1
 40.5, //2
 54.0, //3
 72.0, //4
 85.0, //5
 86.0, //6
 73.5, //7
 54.0, //8
 44.5, //9
 54.0, //10
 45.5, //11
 54.0, //12
 57.0, //13
 81.0, //14
 81.0, //15
 54.0, //16
 42.0, //17
 54.0, //18
 46.5, //19
 54.0, //20
 67.5, //21
 67.5, //22
 56.5, //23
 70.0, //24
 54.0, //25
 109.5, //26
 54.0, //27
 83.5, //28
 67.5, //29
 15, //30
 16, //31
 0, //32
};

extern double *weight_atoms(int nratoms,  const int *atomtypes) {
  //Allocates a new array that contain the weight of an array of atomtypes
  double *atomweights = new double[nratoms];
  for (int n = 0; n < nratoms; n++) {
    int typ = atomtypes[n];
    if (typ < 0 || typ > 32) {
      fprintf(stderr, "Atom weighting is incompatible with non-ATTRACT atom type %d\n", typ);
      exit(1);
    }
    atomweights[n] = aweights[typ];
  }
  return atomweights;
}  