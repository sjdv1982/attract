#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

const int ATOMDENSITYGRIDS = 100;

struct AtomDensityGrid {
  float voxelsize;
  int dimx, dimy, dimz;
  float orix, oriy, oriz;
  double *density;
  bool *maskgrid;
  float clash_threshold;
  float forceconstant;
};
AtomDensityGrid grids[ATOMDENSITYGRIDS];
int ngrids = 0;

typedef void (*Applyfunc) (AtomDensityGrid &g, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight);

inline void trilin(AtomDensityGrid &g, Applyfunc apply, double ax, double ay, double az, double weight, double &overlap, double &gradx, double &grady, double &gradz) {
  double wx1=0,wy1=0,wz1=0;
    
  int px0 = floor(ax);  
  int px1 = ceil(ax);
  if (px1 < 0 || px0 >= g.dimx) return;
  if (px0 < 0) wx1 = 1;
  else if (px1 >= g.dimx) wx1 = 0;
  else wx1 = ax - px0;

  int py0 = floor(ay);  
  int py1 = ceil(ay);
  if (py1 < 0 || py0 >= g.dimy) return;
  if (py0 < 0) wy1 = 1;
  else if (py1 >= g.dimy) wy1 = 0;
  else wy1 = ay - py0;

  int pz0 = floor(az);  
  int pz1 = ceil(az);
  if (pz1 < 0 || pz0 >= g.dimz) return;
  if (pz0 < 0) wz1 = 1;
  else if (pz1 >= g.dimz) wz1 = 0;
  else wz1 = az - pz0;
  
  double wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;
  double dimxy = g.dimx*g.dimy;
  double cwx, cwxy, cwxyz;
  int indx,indxy,indxyz; 
  if (wx0 > 0) {
    cwx = wx0;
    indx = px0;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + g.dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy0*wz0,-wx0*wz0,-wx0*wy0, weight);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy0*wz1,-wx0*wz1,wx0*wy0, weight);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + g.dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy1*wz0,wx0*wz0,-wx0*wy1, weight);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,-wy1*wz1,wx0*wz1,wx0*wy1, weight);
      }
    }
  }
  if (wx1 > 0) {
    cwx = wx1;
    indx = px1;
    if (wy0 > 0){
      cwxy = cwx * wy0;
      indxy = indx + g.dimx*py0;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,wy0*wz0,-wx1*wz0,-wx1*wy0, weight);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,wy0*wz0,-wx1*wz1,-wx1*wy0, weight);
      }
    }
    if (wy1 > 0) {
      cwxy = cwx * wy1;
      indxy = indx + g.dimx*py1;      
      if (wz0 > 0) {
        cwxyz = cwxy * wz0;
        indxyz = indxy + dimxy*pz0;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,wy1*wz0,wx1*wz0,-wx1*wy1, weight);
      }
      if (wz1 > 0) {
        cwxyz = cwxy * wz1;
        indxyz = indxy + dimxy*pz1;      
        apply(g,indxyz,cwxyz,overlap,gradx,grady,gradz,wy1*wz1,wx1*wz1,wx1*wy1, weight);
      }
    }
  }
}

void apply_weight(AtomDensityGrid &g, int indxyz, double cwxyz, double &dummy, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight) {
  g.density[indxyz] += weight * cwxyz;         
}

void apply_clash(AtomDensityGrid &g, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight) {
  double w = g.density[indxyz];
  double clash = w/g.clash_threshold-1;
  if (clash <= 0) return; 
  overlap += clash*clash;
  double cgrad = 2*clash*cwxyz;
  gradx -= sx * cgrad;
  grady -= sy * cgrad;
  gradz -= sz * cgrad;
}

void apply_clash_mask(AtomDensityGrid &g, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight) {
  double w = g.density[indxyz];
  double clash = w/g.clash_threshold-g.maskgrid[indxyz];
  if (clash <= 0) return; 
  overlap += clash*clash;
  double cgrad = 2*clash*cwxyz;
  gradx -= sx * cgrad;
  grady -= sy * cgrad;
  gradz -= sz * cgrad;
}

void apply_clash_noforces(AtomDensityGrid &g, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight) {
  double w = g.density[indxyz];
  double clash = w/g.clash_threshold-1;
  if (clash <= 0) return; 
  overlap += clash*clash;
}

void apply_clash_mask_noforces(AtomDensityGrid &g, int indxyz, double cwxyz, double &overlap, double &gradx, double &grady, double &gradz, double sx, double sy, double sz, double weight) {
  double w = g.density[indxyz];
  double clash = w/g.clash_threshold-g.maskgrid[indxyz];
  if (clash <= 0) return; 
  overlap += clash*clash;
}

static void error(const char *filename) {
  fprintf(stderr, "Reading error in mask file %s\n", filename);
  exit(1);
}

extern "C" void define_atomdensitymask_(char *maskfile, float forceconstant) {
  if (ngrids == ATOMDENSITYGRIDS) {
    fprintf(stderr, "Too many atom density grids\n");
    exit(1);
  }
  FILE *f = fopen(maskfile, "rb");
  if (!f) {
    fprintf(stderr, "Mask file %s does not exist\n", maskfile);
    exit(1);
  }
  AtomDensityGrid &g = grids[ngrids];
  char attractmask[12];
  attractmask[11]=0;
  int read;
  read = fread(attractmask, 11, 1, f);
  if (!read) error(maskfile);
  if (strcmp(attractmask, "ATTRACTMASK")) {
    fprintf(stderr, "File %s is not a mask file\n", maskfile);
    exit(1);    
  }
  read = fread(&g.voxelsize, sizeof(float), 1, f);
  if (!read) error(maskfile);
  read = fread(&g.orix, sizeof(float), 1, f);
  if (!read) error(maskfile);  
  read = fread(&g.oriy, sizeof(float), 1, f);
  if (!read) error(maskfile);  
  read = fread(&g.oriz, sizeof(float), 1, f);
  if (!read) error(maskfile);      
  read = fread(&g.dimx, sizeof(int), 1, f);
  if (!read) error(maskfile);
  read = fread(&g.dimy, sizeof(int), 1, f);
  if (!read) error(maskfile);
  read = fread(&g.dimz, sizeof(int), 1, f);
  if (!read) error(maskfile);
  g.clash_threshold = 0.990 * g.voxelsize*g.voxelsize*g.voxelsize; //average protein packing density is 0.824 Da/A**3 => 0.990 Da/A**3 with 20 % margin
  int nvox = g.dimx*g.dimy*g.dimz;
  g.density = new double[nvox];  
  g.maskgrid = new bool[nvox];
  if ((!g.density) || (!g.maskgrid)) {
    fprintf(stderr, "Cannot allocate memory for atom density grid of dimension %d x %d x %d\n", g.dimx,g.dimy,g.dimz);
    exit(1);    
  }
  read = fread(g.maskgrid, nvox*sizeof(bool), 1, f);
  if (!read) error(maskfile);
  g.forceconstant = forceconstant;
  ngrids++;
}

extern "C" void define_atomdensitygrid_(float voxelsize, int dimension, float forceconstant) {
  if (ngrids == ATOMDENSITYGRIDS) {
    fprintf(stderr, "Too many atom density grids\n");
    exit(1);
  }
  AtomDensityGrid &g = grids[ngrids];
  g.voxelsize = voxelsize;
  g.dimx = dimension;
  g.dimy = dimension;
  g.dimz = dimension;
  g.orix = -0.5 * dimension * g.voxelsize;
  g.oriy = -0.5 * dimension * g.voxelsize;
  g.oriz = -0.5 * dimension * g.voxelsize;  
  g.clash_threshold = 0.990 * voxelsize*voxelsize*voxelsize; //average protein packing density is 0.824 Da/A**3 => 0.990 Da/A**3 with 20 % margin
  g.density = new double[dimension*dimension*dimension];  
  if (!g.density) {
    fprintf(stderr, "Cannot allocate memory for atom density grid of dimension %d\n", dimension);
    exit(1);    
  }
  g.maskgrid = NULL;
  g.forceconstant = forceconstant;
  ngrids++;
}

extern double *weight_atoms(int nratoms,  const int *atomtypes);

extern "C" void atomdensitygridenergy_(double &energy, const int &nratoms, const double *atoms, const double *atoms0, const int *atomtypes, const int &nlig, const int *natomlig, double *forces, const int &update_forces) {  
  energy = 0;
  Applyfunc app;    
  double *atomweights = weight_atoms(nratoms, atomtypes);
  for (int i=0;i < ngrids;i++) {
    AtomDensityGrid &g = grids[i];
    int gridsize = g.dimx*g.dimy*g.dimz;
    memset(g.density, 0, gridsize * sizeof(double));    
    
    //Calculate all the densities
    for (int n=0;n<nratoms;n++) {
      double ax = (atoms[3*n]-g.orix)/g.voxelsize;    
      double ay = (atoms[3*n+1]-g.oriy)/g.voxelsize;
      double az = (atoms[3*n+2]-g.oriz)/g.voxelsize;
      //printf("%.3f %.3f %.3f\n", ax, ay, ax);
      
      app = apply_weight;
      double atomweight = atomweights[n];
      double dmmy;
      trilin(g,app,ax,ay,az, atomweight, dmmy, dmmy, dmmy, dmmy);         
    }   
    
    //Detect clashes (regions with too high atom densities), 
    // and apply repulsive forces to atoms in those regions
    for (int n=0;n<nratoms;n++) {
      double ax = (atoms[3*n]-g.orix)/g.voxelsize;    
      double ay = (atoms[3*n+1]-g.oriy)/g.voxelsize;
      double az = (atoms[3*n+2]-g.oriz)/g.voxelsize;
      double totclash = 0;
      double gradx = 0, grady = 0, gradz = 0;
      
      if (g.maskgrid) {
        app = apply_clash_mask;
        if (!update_forces) app = apply_clash_mask_noforces;        
      }  
      else {
        app = apply_clash;
        if (!update_forces) app = apply_clash_noforces;                
      }  
      trilin(g,app, ax,ay,az, 0, totclash, gradx,grady,gradz);

      energy += totclash * g.forceconstant;
      if (update_forces) {
        gradx *= g.forceconstant;
        grady *= g.forceconstant;
        gradz *= g.forceconstant;
        forces[3*n] += gradx;
        forces[3*n+1] += grady;
        forces[3*n+2] += gradz;
      }
    }
  }
  delete[] atomweights;
}
