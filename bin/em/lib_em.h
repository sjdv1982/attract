#include "stdlib.h"
#include "math.h"

#define g1idz(k, j, i) g_extx*g_exty*k+g_extx*j+i
typedef double Coor[3];

#ifdef __cplusplus
extern "C" {
#endif

double gridify (
 Coor *pdb2, int num2, 
 double *phi2,  int nvox2, 
 double width, int g_extx, int g_exty, int extz,
 double minx, double miny, double minz
);

double *calc_kernel(double width, double reso, double kampl, unsigned int *pg_extx);

void apply_kernel(
  /*input map*/ const double *phi, unsigned int nvox, unsigned int g_extx, unsigned int g_exty, unsigned int extz,
  /*kernel map*/ const double *phi2, unsigned int nvox2, unsigned int g_extx2, unsigned int g_exty2, unsigned int extz2,
  /*output map*/ double *phi3
);

#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif
