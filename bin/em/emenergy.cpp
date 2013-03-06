#include "em.h"

Map maps[100];
int nrmaps = 0;

typedef void (*Applyfunc) (Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz);

void apply_overlap(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  overlap += cwxyz * m.precomp_overlap[indxyz];
  gradx += cwxyz * m.precomp_grad[3*indxyz];
  grady += cwxyz * m.precomp_grad[3*indxyz+1];
  gradz += cwxyz * m.precomp_grad[3*indxyz+2];
}

void apply_overlap_noforces(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  overlap += cwxyz * m.precomp_overlap[indxyz];
}

void apply_overlap_and_weight(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  overlap += cwxyz * m.precomp_overlap[indxyz];
  gradx += cwxyz * m.precomp_grad[3*indxyz];
  grady += cwxyz * m.precomp_grad[3*indxyz+1];
  gradz += cwxyz * m.precomp_grad[3*indxyz+2];
  m.weights[indxyz] += transw(cwxyz);          
}

void apply_overlap_and_weight_noforces(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  overlap += cwxyz * m.precomp_overlap[indxyz];
  m.weights[indxyz] += transw(cwxyz);          
}

void apply_clash(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  double w = m.weights[indxyz];
  double clash = w - m.clash_threshold;
  if (clash <= 0) return; 
  overlap += clash*clash;
  double cgrad = 2*clash*gradw(cwxyz);
  gradx -= sx * cgrad;
  grady -= sy * cgrad;
  gradz -= sz * cgrad;
}

void apply_clash_noforces(Map &m, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz) {
  double w = m.weights[indxyz];
  double clash = w - m.clash_threshold;
  if (clash <= 0) return; 
  overlap += clash*clash;
}

inline void trilin(Map &m, Applyfunc apply, double ax, double ay, double az, double &overlap, double &gradx, double &grady, double &gradz) {
  double wx1=0,wy1=0,wz1=0;
  
  
  int px0 = floor(ax);  
  int px1 = ceil(ax);
  if (px1 < 0 || px0 >= m.dimx) return;
  if (px0 < 0) wx1 = 1;
  else if (px1 >= m.dimx) wx1 = 0;
  else wx1 = ax - px0;

  int py0 = floor(ay);  
  int py1 = ceil(ay);
  if (py1 < 0 || py0 >= m.dimy) return;
  if (py0 < 0) wy1 = 1;
  else if (py1 >= m.dimy) wy1 = 0;
  else wy1 = ay - py0;

  int pz0 = floor(az);  
  int pz1 = ceil(az);
  if (pz1 < 0 || pz0 >= m.dimz) return;
  if (pz0 < 0) wz1 = 1;
  else if (pz1 >= m.dimz) wz1 = 0;
  else wz1 = az - pz0;
  
  double wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;
  double dimxy = m.dimx*m.dimy;
  double cwx, cwxy, cwxyz;
  int indx,indxy,indxyz; 
  if (wx0 > 0) {
    cwx = wx0;
    indx = px0;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + m.dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy0*wz0,-wx0*wz0,-wx0*wy0);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy0*wz1,-wx0*wz1,wx0*wy0);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + m.dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy1*wz0,wx0*wz0,-wx0*wy1);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy1*wz1,wx0*wz1,wx0*wy1);
      }
    }
  }
  if (wx1 > 0) {
    cwx = wx1;
    indx = px1;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + m.dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,wy0*wz0,-wx1*wz0,-wx1*wy0);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,wy0*wz0,-wx1*wz1,-wx1*wy0);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + m.dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,wy1*wz0,wx1*wz0,-wx1*wy1);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(m,indxyz,cwxyz,overlap,gradx,grady,gradz,wy1*wz1,wx1*wz1,wx1*wy1);
      }
    }
  }
}


extern "C" void read_densitymaps_(char *densitymap0, float resolution, int len_densitymap) {
  char *densitymap = new char[len_densitymap+1];
  densitymap[len_densitymap] = 0;
  memcpy(densitymap, densitymap0, len_densitymap);
  if (nrmaps == 0) memset(maps, 0, 100 * sizeof(Map));
  {
  Map &m = maps[nrmaps+1];
  m.mode = EM_MAPMODE_CLASH;
  m.clash_threshold = 16 * m.situs_width * m.situs_width * m.situs_width / 1000;
  m.clash_weight = 0.1;
  read_vol(densitymap, &m.situs_width, &m.situs_origx,&m.situs_origy,&m.situs_origz,&m.dimx,&m.dimy,&m.dimz,&m.densities);          
  int gridsize = m.dimx * m.dimy * m.dimz;  
  m.weights = new double[gridsize];  
  }
  
  {
  Map &m = maps[nrmaps];
  memcpy(&m,  &maps[nrmaps+1], sizeof(Map)); //copy the clash map
  m.mode = EM_MAPMODE_OVERLAP_AND_WEIGHT;
  m.resolution = resolution * 8;
  m.sigma = 0.5 * m.resolution / sqrt(3);  
  m.nrmaxatoms = 0;
  m.emweight = 10000000;
  m.overlapmargin = 1.0;
  m.fac = 4.0;
  m.base = exp(m.fac);
  precompute(m, 1); //HAS to be 1!
  }  
  nrmaps+=2;
  //2  
  {
  Map &m = maps[nrmaps]; 
  memcpy(&m,  &maps[nrmaps-2], sizeof(Map)); //copy the overlap+weight map
  m.resolution = resolution * 4;
  m.sigma = 0.5 * m.resolution / sqrt(3); 
  m.emweight = 1000000;
  m.overlapmargin = 1.0;
  precompute(m, 1);
  nrmaps++;
  }

  //3  
  {
  Map &m = maps[nrmaps]; 
  memcpy(&m,  &maps[nrmaps-1], sizeof(Map));  //copy the previous overlap map
  m.resolution = resolution * 2.5;
  m.sigma = 0.5 * m.resolution / sqrt(3); 
  m.emweight = 100000;
  m.overlapmargin = 0.85;
  precompute(m, 1);
  nrmaps++;
  }
  
  //4
  {
  Map &m = maps[nrmaps]; 
  memcpy(&m,  &maps[nrmaps-1], sizeof(Map));  //copy the previous overlap map
  m.resolution = resolution * 2;
  m.sigma = 0.5 * m.resolution / sqrt(3); 
  m.emweight = 100000;
  m.overlapmargin = 0.85;
  precompute(m, 1);
  nrmaps++;
  }

  //5
  {
  Map &m = maps[nrmaps]; 
  memcpy(&m,  &maps[nrmaps-1], sizeof(Map));  //copy the previous overlap map
  m.resolution = resolution * 1.5;
  m.sigma = 0.5 * m.resolution / sqrt(3); 
  m.emweight = 100000;  
  m.overlapmargin = 0.85;
  precompute(m, 1);
  nrmaps++;
  }
  
  delete[] densitymap;
}

double emenergy (Map &m, int nratoms, const double *atoms, const double *atoms0, const int *atomtypes, int nlig, const int *natomlig, double *forces, bool update_forces, double &totgradx, double &totgrady, double &totgradz) {  
  double energy = 0;
  int n;    
  Applyfunc app;    
  if (m.mode == EM_MAPMODE_OVERLAP_AND_WEIGHT || m.mode == EM_MAPMODE_OVERLAP)  
  {
    if (m.emweight == 0) return 0;
    
    if (m.nrmaxatoms == 0) {
      //TODO?: supply atomtypes and properly weight the atoms
      em_maxoverlap(m, nratoms, atoms0, nlig, natomlig);
    }
    totgradx = 0; totgrady = 0; totgradz = 0;
    if (m.mode == EM_MAPMODE_OVERLAP_AND_WEIGHT) {
      int gridsize = m.dimx*m.dimy*m.dimz;
      memset(m.weights, 0, gridsize * sizeof(double));
    }      
    for (n=0;n<nratoms;n++) {
      double ax = (atoms[3*n] - m.situs_origx)/m.situs_width;    
      double ay = (atoms[3*n+1] - m.situs_origy)/m.situs_width;
      double az = (atoms[3*n+2] - m.situs_origz)/m.situs_width;
      double totoverlap = 0;
      double gradx = 0, grady = 0, gradz = 0;
      
      /*
      ///
      double voxels_per_sigma = m.sigma / m.situs_width;
      calc_emdensity(voxels_per_sigma, m.base,m.fac,
       m.dimx, m.dimy, m.dimz, m.densities,
       ax,ay,az, 
       totoverlap, gradx, grady, gradz);
      goto x;
      ///
      */
      
      
      if (m.mode == EM_MAPMODE_OVERLAP) {
        app = apply_overlap;
        if (!update_forces) app = apply_overlap_noforces;
      }  
      else {
        app = apply_overlap_and_weight;
        if (!update_forces) app = apply_overlap_and_weight_noforces;
      }  
      trilin(m,app,ax,ay,az,totoverlap, gradx,grady,gradz);
//x:
      double max_overlap = m.maxatoms[n];  
      
      double gradscale = 2*m.emweight/(-max_overlap*m.situs_width);
      double lackfac = max_overlap/m.maxmax;
      double totoverlapf = totoverlap / max_overlap;
      //printf("OVER %d %.3f %.3f %.3f\n", m.mode, max_overlap, totoverlap, totoverlapf);
      double curr_energy = 0;
      if (totoverlapf < m.overlapmargin) {
        double lack = (m.overlapmargin - totoverlapf) * lackfac;
        curr_energy = m.emweight * lack*lack;
        energy += curr_energy;
        double currgradscale = lack * gradscale;
        if (update_forces == 1) {
          gradx *= currgradscale;
          grady *= currgradscale;
          gradz *= currgradscale;
          totgradx += gradx; totgrady += grady; totgradz += gradz;        
          forces[3*n] += gradx;
          forces[3*n+1] += grady;
          forces[3*n+2] += gradz;
        }
      }
    }
  }  
  else if (m.mode == EM_MAPMODE_CLASH) { 
    for (n=0;n<nratoms;n++) {
      double ax = (atoms[3*n] - m.situs_origx)/m.situs_width;    
      double ay = (atoms[3*n+1] - m.situs_origy)/m.situs_width;
      double az = (atoms[3*n+2] - m.situs_origz)/m.situs_width;
      double totclash = 0;
      double gradx = 0, grady = 0, gradz = 0;

      app = apply_clash;
      if (!update_forces) app = apply_clash_noforces;
      trilin(m,app, ax,ay,az,totclash, gradx,grady,gradz);

      energy += totclash * m.clash_weight;
      if (update_forces) {
        gradx *= m.clash_weight;
        grady *= m.clash_weight;
        gradz *= m.clash_weight;
        totgradx += gradx; totgrady += grady; totgradz += gradz;        
        forces[3*n] += gradx;
        forces[3*n+1] += grady;
        forces[3*n+2] += gradz;
      }
    }
  }
  return energy;  
}


extern "C" void emenergy_(double &energy, const int &nratoms, const double *atoms, const double *atoms0, const int *atomtypes, const int &nlig, const int *natomlig, double *forces, const int &update_forces) {  

  energy = 0;
  for (int m = 0; m < nrmaps; m++) {
    double dumx,dumy,dumz;
    double denergy = emenergy(maps[m], nratoms, atoms, atoms0, atomtypes, nlig, natomlig, forces, update_forces,dumx,dumy,dumz);
    energy += denergy;
  }
}
