#include "em.h"
  /*
  const int maxmaxatoms = 100000;
  double maxatoms[maxmaxatoms];
  int nrmaxatoms = 0;
  double maxmax = 0;
  //...
  m.maxatoms = new double[nrmaxatoms];
  memcpy(m.maxatoms, maxatoms, nrmaxatoms*sizeof(double));
  m.nrmaxatoms = nrmaxatoms;
  m.maxmax = maxmax;  
  */

void calc_overlap(Map &m, double dissq, double dx, double dy, double dz, double &totoverlap, int index, double &gradx, double &grady, double &gradz) {
  double density = m.densities[index];
  double overlap = density*pow(m.base, -dissq);
  totoverlap += overlap;
  double fgrad = m.fac * -2 * overlap;
  gradx += dx*fgrad;
  grady += dy*fgrad;
  gradz += dz*fgrad; 
}

double em_maxoverlap_ligand(Map &m, int nratoms, const double *atoms, double *maxoverlaps) {
  double maxdis = sigma_threshold;
  double maxdissq = maxdis * maxdis;
  double mindis = m.situs_width / m.sigma; 
  double mindissq = mindis * mindis;
  double maxmax = 0;
  for(int i = 0; i < nratoms; i++) {
    double &maxoverlap = maxoverlaps[i];
    const double *a1 = &atoms[3*i];
    for(int j = 0; j < nratoms; j++) {
      const double *a2 = &atoms[3*j];
      double dissq = 0;
      for (int k = 0; k < 3; k++) {
        double dk = a2[k] - a1[k];
        dissq += dk *dk;
      }
      if (dissq > maxdissq) continue;
      if (dissq < mindissq) dissq = mindissq;
      double overlap = pow(m.base, -dissq);
      maxoverlap += overlap; 
    } 
    //printf("MAXOVERLAP %.5f\n", maxoverlap);
    if (maxoverlap > maxmax) maxmax = maxoverlap;
  }
  return maxmax;
}
void em_maxoverlap(Map &m, int nratoms, const double *atoms, int nlig, const int *natomlig) {
  m.maxatoms = new double[nratoms];
  memset(m.maxatoms, 0, nratoms*sizeof(double));
  int pos = 0;
  double *atoms2 = new double[3*nratoms];
  int n;
  //printf("%.3f\n", m.sigma);
  for (n = 0; n < 3*nratoms; n++) { 
    atoms2[n] = atoms[n] / m.sigma;
  }
  double maxmax = 0;
  for (n = 0; n < nlig; n++) { 
    double lmaxmax = em_maxoverlap_ligand(m, natomlig[n], &atoms2[3*pos], &m.maxatoms[pos]);
    if (lmaxmax > maxmax) maxmax = lmaxmax;
    pos += natomlig[n];
  }
  m.maxmax = maxmax;
  m.nrmaxatoms = nratoms;
  delete[] atoms2;
}
  
