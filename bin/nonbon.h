#include "grid.h"
#include <cmath>

inline void nonbon(int iab, double welwel, double rc, double ac, double emin, double rmin2, int ivor, double dsq, double rr2, double dx, double dy, double dz, int potshape, double fswi,
double &energy, Coor &grad) 

{
  double alen = welwel * ac; 
  double rlen = welwel * rc; 
    
  double rr23 = rr2 * rr2 * rr2;
  double rrd, shapedelta;
    
  if (potshape==8) {
    rrd = rr2;
    shapedelta = 2;
  }
  else if (potshape==12) {
    rrd = rr23;
    shapedelta = 6;    
  }
  double rep = rlen * rrd;
  double vlj = (rep-alen)*rr23; 
  
  if (dsq < rmin2) {
    energy = vlj + fswi * (ivor-1) * emin;
    if (iab) {
      double fb=fswi*6.0*vlj+shapedelta*(rep*rr23);
      grad[0] = fb * dx;
      grad[1] = fb * dy;
      grad[2] = fb * dz;
    }
  }
  else {
    energy = fswi*ivor * vlj;
    if (iab) {
      double fb=fswi*6.0*vlj+shapedelta*(rep*rr23);
      grad[0] = ivor * fb * dx;
      grad[1] = ivor * fb * dy;
      grad[2] = ivor * fb * dz;
    }
  }
  
}

//TODO: update with swi and potshape!
inline void nonbon_nograd(double welwel, double rc, double ac, double emin, double rmin2, int ivor, double dsq, double rr2, 
double &energy) 

{
  double alen = welwel * ac; 
  double rlen = welwel * rc; 
    
  double rr23 = rr2 * rr2 * rr2;
  double rep = rlen * rr2;
  double vlj = (rep-alen)*rr23; 
  
  if (dsq < rmin2) {
    energy = vlj + (ivor-1) * emin;
  }
  else {
    energy = ivor * vlj;
  }
  
}

const double felec = double(332.053986)/15.0;
const double felecsqrt = sqrt(felec);

inline void elec(int iab, double charge, 
double rr2, double dx, double dy, double dz, 
double &energy, Coor &grad) {

  double et = charge * rr2;
  energy = et;
  if (iab) {
    grad[0] = 2 * et * dx;
    grad[1] = 2 * et * dy;
    grad[2] = 2 * et * dz;
  }
}

inline void elec_nograd(double charge, 
double rr2, double &energy) {

  double et = charge * rr2;
  energy = et;
}

