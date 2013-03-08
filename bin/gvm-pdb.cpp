//Calculates Gradient Vector Matching score on a PDB file

//usage: ./gvm-pdb map.vol <gradient threshold> <PDB file>

#include "max.h"
#include "em/lib_em.h"
#include "em/lib_corr.h"

/*support for non-reduced PDBs*/
#include "state.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

extern void read_pdb2(
  FILE *fil, Coor *&x, 
  char **&pdbstrings, bool *&pdblayout,
  int &coorcounter, int &linecounter
);

const double p = -1.0f; //sign swaps, because we are *convoluting* with a kernel; an atom at X=10 means a negative X gradient at X=11!
const double gvm_kernel_x[] = {
 -p, 0, p,
 -p, 0, p,
 -p, 0, p,

 -p, 0, p,
 -p, 0, p,
 -p, 0, p,

 -p, 0, p,
 -p, 0, p,
 -p, 0, p
};
const double gvm_kernel_y[] = {
 -p,-p,-p,
  0, 0, 0,
  p, p, p,

 -p,-p,-p,
  0, 0, 0,
  p, p, p,

 -p,-p,-p,
  0, 0, 0,
  p, p, p
};
const double gvm_kernel_z[] = {
 -p,-p,-p,
 -p,-p,-p,
 -p,-p,-p,

  0, 0, 0,
  0, 0, 0,
  0, 0, 0,

  p, p, p,
  p, p, p,
  p, p, p  
};

extern "C" void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi);

extern int cartstate_new(int argc, char *argv[],bool single=0);

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[3][MAXLIG]);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_f_write_pdb_(
  const int &handle,
  int &nlig, int *&kai, char4 *&tyi, char4 *&rgi, int *&iei, double *&x,
  int *&iaci, double *&xlai, int *&icop, double *&we, int *&ieins);

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *(&nhm), int *(&ieins), double *(&eig), double *(&pivot), double *(&xb), double *(&x),double *(&xori), double *(&xori0));

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nlm));
extern "C" void cartstate_get_pivot_(const int &handle,double *&pivot);
extern "C" void cartstate_get_nrens_(const int &handle,int *&nrens);

extern "C" void cartstate_get_ensd_(const int &handle,
  const int &ligand,
  const int &ens,
  double *&ensd,
  const double &morph,
  double &cmorph,
  double *&cmorphd
  );
  
extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, coors2 &locrests, dof2 &morph, modes2 &dlig, const int &nlig, const int *nhm, const int *nrens0, const int *morphing, const int *has_locrests, int &seed, char *&label, const int &all_labels, int f_len);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

extern "C" void rotate_(const int &maxlig,const int &max3atom,double (&rotmat)[9],const double &xa,const double &ya,const double &za,
double *pivot,
int &ijk,int *ieins, double *x);

extern "C" void deform_(const int &maxlig,const int &max3atom, 
const int &totmax3atom, const int &maxatom,const int &maxmode,
int &ens, double *ensdp, const double &cmorph, const double *cmorphdp, 
double (&dligp)[MAXMODE], 
int *nhm,int &ijk,int *ieins,double *eig,double *xb,double *x,double *xori,double *xori0, const int &do_morph); 

extern void read_ens(int cartstatehandle, int ligand, char *ensfile, bool strict, bool morphing);

CartState &cartstate_get(int handle);

void usage() {
  fprintf(stderr, "usage: gvm-pdb map.vol <gradient threshold> <PDB file>\n");
  exit(1);
}

bool exists(const char *f) {
  FILE *fil = fopen(f, "r");
  if (fil == NULL) return 0;
  else {
    fclose(fil);
    return 1;
  }
}

int main(int argc, char *argv[]) {
  int i;
  if (argc !=  4) {
    fprintf(stderr, "Wrong number of arguments\n"); usage();
  }

  for (i = 1; i < argc; i++) {
    if (i == 2) continue;
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
  }

  char **pdbstrings[MAXLIG]; bool *pdblayout[MAXLIG]; int linecounter[MAXLIG];
  Coor *pdb;
  int nratoms;
  FILE *fil = fopen(argv[3], "r");
  read_pdb2(fil,pdb,pdbstrings[i],pdblayout[i],nratoms,linecounter[i]);
 
  double width, minx, miny, minz;
  unsigned int g_extx; unsigned int g_exty; unsigned int extz;
  double *mapdata; 

  read_vol(argv[1], &width, &minx, &miny, &minz, &g_extx, &g_exty, &extz, &mapdata);
  unsigned int nvox = g_extx * g_exty * extz;
  double threshold = atof(argv[2]);

  double *mapdata_xyz = (double *) malloc(3*nvox*sizeof(double));     
  double *mapdata_x = mapdata_xyz;
  apply_kernel(
    mapdata, nvox, g_extx, g_exty, extz,
    gvm_kernel_x, 27, 3, 3, 3,
    mapdata_x
  );

  double *mapdata_y = &mapdata_xyz[nvox];
  apply_kernel(
    mapdata, nvox, g_extx, g_exty, extz,
    gvm_kernel_y, 27, 3, 3, 3,
    mapdata_y
  );

  double *mapdata_z = &mapdata_xyz[2*nvox];
  apply_kernel(
    mapdata, nvox, g_extx, g_exty, extz,
    gvm_kernel_z, 27, 3, 3, 3,
    mapdata_z
  );

  double *mappdb = (double *) malloc(nvox*sizeof(double));
  double *mappdb_xyz = (double *) malloc(3*nvox*sizeof(double));     
   
  gridify(
   pdb, nratoms, 
   mappdb, nvox, 
   width, g_extx, g_exty, extz,
   minx, miny, minz  
  );

  double *mappdb_x = mappdb_xyz;
  apply_kernel(
    mappdb, nvox, g_extx, g_exty, extz,
    gvm_kernel_x, 27, 3, 3, 3,
    mappdb_x
  );

  double *mappdb_y = &mappdb_xyz[nvox];
  apply_kernel(
    mappdb, nvox, g_extx, g_exty, extz,
    gvm_kernel_y, 27, 3, 3, 3,
    mappdb_y
  );

  double *mappdb_z = &mappdb_xyz[2*nvox];
  apply_kernel(
    mappdb, nvox, g_extx, g_exty, extz,
    gvm_kernel_z, 27, 3, 3, 3,
    mappdb_z
  );
     
  double sumx = 0, sumy = 0, sumxx = 0, sumxy = 0, sumyy = 0; 
  int count = 0;
  for (int dim = 0; dim < 3; dim++) {
    double *mapd = mapdata_xyz + nvox * dim;
    double *mapw = mappdb_xyz + nvox * dim;
    for (int z = 1; z < extz-1; z++) {
      for (int y = 1; y < g_exty-1; y++) {
        for (int x = 1; x < g_extx-1; x++) {      
          int n = g_exty*g_extx*z + g_extx * y + x;
          double d = mapd[n];
          if (fabs(d) < threshold) continue;
          double w = mapw[n];
          sumx += d; sumxx += d*d;
          sumy += w; sumyy += w*w;
          sumxy += d * w;
          count++;
        }
      }
    }
  }
  double Sxx = sumxx - sumx * sumx / count;
  double Sxy = sumxy - sumx * sumy / count;
  double Syy = sumyy - sumy * sumy / count;
  double r = Sxy/sqrt(Sxx*Syy);

  printf("%.6f\n",r);
}
