#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <memory.h>

const double sigma_threshold = 3;

inline double transw(double weight) {
  if (weight > 0.5) return 0.5-(1-weight)*(1-weight);
  return weight*weight;
}

inline double gradw(double weight) {
  if (weight > 0.5) return 2 * (1-weight);
  return 2*weight;
}

enum MapMode {
  EM_MAPMODE_OVERLAP,
  EM_MAPMODE_CLASH,
  EM_MAPMODE_OVERLAP_AND_WEIGHT,
};

struct Map {
  //General parameters
  unsigned int dimx;
  unsigned int dimy;
  unsigned int dimz;
  double situs_origx, situs_origy, situs_origz, situs_width;
  float resolution;

  MapMode mode;
  
  
  //parameters for overlap mode
  double overlapmargin; //fraction of maximum overlap that must be reached;less will give an energy penalty
  double emweight;      // weight of the energy of this map 
  double *maxatoms;
  int nrmaxatoms;
  double maxmax; //maximum value of maxatoms
  double *precomp_overlap;
  double *precomp_grad;

  //parameters for overlap mode
  double *densities;
  double clash_threshold;
  double clash_weight;

  //temporary values
  double sigma;
  double base;
  double fac;
  double *weights;  //transformed

};

extern "C" void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi);

extern void em_maxoverlap(Map &m, int nratoms, const double *atoms, int nlig, const int *natomlig);
extern void precompute(Map &m, int sampling);

extern void calc_emdensity(double voxels_per_sigma, double base, double fac,
int dimx, int dimy, int dimz,
double *densities,
double ax, double ay, double az, 
double &totoverlap, double &gradx, double &grady, double &gradz);
