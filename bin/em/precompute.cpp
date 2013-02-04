#include "em.h"
#include <cstdio>
#include <cmath>
#include <memory.h>

inline void calc_overlap(double dissq, double base, double fac, double sigma, double dx, double dy, double dz, double &totoverlap, double density, double &gradx, double &grady, double &gradz) {
  double overlap = density*pow(base, -dissq);
  totoverlap += overlap;
  /*
  double dis = sqrt(dissq); //in sigma!
  double grad = fac * -2 * dis * overlap; /gives the gradient per sigma unit!
  double gradang = grad / sigma; //gradient per angstrom unit
  double fgrad = gradang / dis;  //(gradx = dx/dis * gradang; dx/dis is the same in Angstroms or in Sigma)
  */  
  double fgrad = fac * -2 * overlap / sigma;
  gradx += dx*fgrad;
  grady += dy*fgrad;
  gradz += dz*fgrad; 
}

void calc_emdensity(double voxels_per_sigma, double base, double fac, double sigma,
int dimx, int dimy, int dimz,
double *densities,
double ax, double ay, double az, 
double &totoverlap, double &gradx, double &grady, double &gradz) 
{

  double maxdist, maxdistsq;
  maxdist = sigma_threshold;
  maxdistsq = maxdist * maxdist;
  
  double disfac = 1.0/voxels_per_sigma;
  double dsqfac = disfac*disfac;

  int px0 = floor(ax);
  if (px0 < 0) px0=0;  
  if (px0 >= dimx) px0 = dimx-1;
  double dx0 = ax - px0;
  double dx0sq = dx0*dx0;
  int px1 = ceil(ax);
  if (px1 < 0) px1=0;  
  if (px1 == px0) px1++;
  if (px1 >= dimx) px1 = dimx-1;
  double dx1 = px1-ax;
  double dx1sq = dx1*dx1;

  int py0 = floor(ay);
  if (py0 < 0) py0=0;  
  if (py0 >= dimy) py0 = dimy-1;
  double dy0 = ay - py0;
  double dy0sq = dy0*dy0;
  int py1 = ceil(ay);
  if (py1 < 0) py1=0;    
  if (py1 == py0) py1++;    
  if (py1 >= dimy) py1 = dimy-1;
  double dy1 = py1-ay;
  double dy1sq = dy1*dy1;

  int pz0 = floor(az);
  if (pz0 < 0) pz0=0;  
  if (pz0 >= dimz) pz0 = dimz-1;
  double dz0 = az - pz0;
  double dz0sq = dz0*dz0;
  int pz1 = ceil(az);
  if (pz1 < 0) pz1=0;     
  if (pz1 == pz0) pz1++;    
  if (pz1 >= dimz) pz1 = dimz-1;
  double dz1 = pz1-az;
  double dz1sq = dz1*dz1;

  //positive x loop
  double dsqx = dsqfac * dx1sq;
  double dx = dx1;
  int x;
  for (x = px1; x < dimx; x++) {
    if (dsqx > maxdistsq) break;

    //positive y loop
    double dsqxy = dsqx + dsqfac * dy1sq;
    double dy = dy1;
    int y;
    for (y = py1; y < dimy; y++){
      if (dsqxy > maxdistsq) break;

      //positive z loop
      double dsqxyz = dsqxy + dsqfac * dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < dimz; z++){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma,dx,dy,dz,totoverlap, density,gradx,grady,gradz);

        //increment positive z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //negative z loop
      dsqxyz = dsqxy + dsqfac * dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, dx,dy,-dz,totoverlap, density,gradx,grady,gradz);

        //increment negative z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //increment positive y loop
      dsqxy += dsqfac * (2 *dy + 1);
      dy++;        
    }


    //negative y loop
    dsqxy = dsqx + dsqfac * dy0sq;
    dy = dy0;
    for (y = py0; y >= 0; y--){
      if (dsqxy > maxdistsq) break;

      //positive z loop
      double dsqxyz = dsqxy + dsqfac * dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < dimz; z++){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, dx,-dy,dz,totoverlap, density,gradx,grady,gradz);

        //increment positive z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //negative z loop
      dsqxyz = dsqxy + dsqfac * dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, dx,-dy,-dz,totoverlap, density,gradx,grady,gradz);

        //increment negative z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //increment negative y loop
      dsqxy += dsqfac * (2 *dy + 1);
      dy++;        
    }


    //increment positive x loop
    dsqx += dsqfac * (2 *dx + 1);
    dx++;
  }

  //negative x loop
  dsqx = dsqfac * dx0sq;
  dx = dx0;
  for (x = px0; x >= 0; x--) {
    if (dsqx > maxdistsq) break;

    //positive y loop
    double dsqxy = dsqx + dsqfac * dy1sq;
    double dy = dy1;
    int y;
    for (y = py1; y < dimy; y++){
      if (dsqxy > maxdistsq) break;

      //positive z loop
      double dsqxyz = dsqxy + dsqfac * dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < dimz; z++){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, -dx,dy,dz,totoverlap, density,gradx,grady,gradz);

        //increment positive z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //negative z loop
      dsqxyz = dsqxy + dsqfac * dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, -dx,dy,-dz,totoverlap, density,gradx,grady,gradz);

        //increment negative z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //increment positive y loop
      dsqxy += dsqfac * (2 *dy + 1);
      dy++;        
    }


    //negative y loop
    dsqxy = dsqx + dsqfac * dy0sq;
    dy = dy0;
    for (y = py0; y >= 0; y--){
      if (dsqxy > maxdistsq) break;

      //positive z loop
      double dsqxyz = dsqxy + dsqfac * dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < dimz; z++){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, -dx,-dy,dz,totoverlap, density,gradx,grady,gradz);

        //increment positive z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //negative z loop
      dsqxyz = dsqxy + dsqfac * dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > maxdistsq) break;
        //dsqxyz is the actual distance squared

        double density = densities[dimx*dimy*z+dimx*y+x];
        calc_overlap(dsqxyz, base,fac,sigma, -dx,-dy,-dz,totoverlap, density,gradx,grady,gradz);

        //increment negative z loop
        dsqxyz += dsqfac * (2 * dz + 1);
        dz++;
      }

      //increment negative y loop
      dsqxy += dsqfac * (2 *dy + 1);
      dy++;        
    }


    //increment negative x loop
    dsqx += dsqfac * (2 *dx + 1);
    dx++;
  }

}

void precompute(Map &m, int sampling) {
  int dimx = m.dimx * sampling;
  int dimy = m.dimy * sampling;
  int dimz = m.dimz * sampling;
  int gridsize = dimx * dimy * dimz;
  m.precomp_overlap = new double[gridsize];
  memset(m.precomp_overlap, 0, gridsize *sizeof(double));
  m.precomp_grad = new double[3*gridsize];
  memset(m.precomp_grad, 0, 3 * gridsize *sizeof(double));
  double voxels_per_sigma = m.sigma / m.situs_width;
  for (int x = 0; x < dimx; x++) {
    double xx = double(x)/sampling;
    for (int y = 0; y < dimy; y++) {
      double yy = double(y)/sampling;
      for (int z = 0; z < dimz; z++) {
        double zz = double(z)/sampling;

        int pos = dimx*dimy*z+dimx*y+x;
        double &overlap = m.precomp_overlap[pos];
        double *grad = &m.precomp_grad[3*pos];
                
        calc_emdensity(voxels_per_sigma, m.base,m.fac,m.sigma,
         m.dimx, m.dimy, m.dimz, m.densities,
         xx,yy,zz, 
         overlap, grad[0], grad[1], grad[2]);
        //printf("PRECOMPUTE %.3f %.3f %.3f %.3f\n", overlap, grad[0], grad[1], grad[2]); 
      }
    }
  }
  m.situs_width /= sampling;
  m.dimx = dimx;
  m.dimy = dimy;
  m.dimz = dimz;     
}
