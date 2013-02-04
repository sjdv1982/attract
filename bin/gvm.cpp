//Calculates Gradient Vector Matching score on an ATTRACT .DAT file

//usage: ./gvm map.vol <gradient threshold> structures.dat receptor.pdb [ligand.pdb] [...] [...] [--modes <modefile>] [--ens/--morph <ligand nr> <ensemble file>]
//  if no ligand.pdb, receptor.pdb is a multi-ligand PDB file 

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

/* DOFs */
static int ens[MAXLIG];
static double morph[MAXLIG];
static double phi[MAXLIG];
static double ssi[MAXLIG];
static double rot[MAXLIG];
static double xa[MAXLIG];
static double ya[MAXLIG];
static double za[MAXLIG];
static double dlig[MAXLIG][MAXMODE];
static int seed;
static char *label;

void usage() {
  fprintf(stderr, "usage: gvm map.vol <gradient threshold> structures.dat receptor.pdb [...] [...] [--modes <modefile>] [--ens/--morph <ligand nr> <ensemble file>]");
  exit(1);
}

extern bool exists(const char *f);

/*
inline void print_contour(int x, int y, int z, double *mapd, double *mapw, int n, double threshold) {
  double d = mapd[n];
  if (fabs(d) < threshold) return;
  double w = mapw[n];
  printf("CONTOUR %d %d %d %.3f %.3f %.3f\n", x,y,z,d-w, d, w);
}
*/

int main(int argc, char *argv[]) {
  int i;
  if (argc < 5) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }

  coors2 locrests;
  int has_locrests[MAXLIG];
  memset(has_locrests, 0, MAXLIG*sizeof(int));
 
  char *modefile = NULL;
  int enscount = 0;
  int ens_ligands[MAXLIG];
  char *ens_files[MAXLIG];
  int morphing[MAXLIG];
  memset(morphing,0,MAXLIG*sizeof(int));
  for (int n = 1; n < argc-1; n++) {
    if (!strcmp(argv[n],"--modes")) {
      modefile = argv[n+1];
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+2 < argc) memcpy(argv2+n,argv+n+2,(argc-n-2)*sizeof(char*));
      argv = argv2;
      argc -= 2;
      break;
    }
  }
  for (int n = 1; n < argc-2; n++) {
    bool has_ens = 0;        
    if (!strcmp(argv[n],"--ens")) {
      has_ens = 1;
    }
    if (!strcmp(argv[n],"--morph")) {
      morphing[enscount] = 1;
      has_ens = 1;
    }    
    if (has_ens) {
      ens_ligands[enscount] = atoi(argv[n+1]);
      ens_files[enscount] = argv[n+2];
      enscount++;
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+3 < argc) memcpy(argv2+n,argv+n+3,(argc-n-3)*sizeof(char*));
      argv = argv2;
      argc -= 3;
      n -= 1;
    }    
  }  
  for (int n = 1; n < argc-1; n++) {
    if (!strcmp(argv[n],"--locrest")) {
      int lig = atoi(argv[n+1]);
      if (lig <= 0 || lig > MAXLIG) {
        fprintf(stderr,"Ligand %d must be larger than 0\n", lig);
        usage();
      }
      has_locrests[lig-1] = 1;
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+2 < argc) memcpy(argv2+n,argv+n+2,(argc-n-2)*sizeof(char*));
      argv = argv2;
      argc -= 2;
      n -= 1;
    }
  }

  for (i = 1; i < argc; i++) {
    if (i == 2) continue;
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
  }

  char **pdbstrings[MAXLIG]; bool *pdblayout[MAXLIG]; int linecounter[MAXLIG];
  //load the Cartesian parameters and get a handle to it
  int cartstatehandle;  
  if (argc == 5) { //one PDB, reduced or non-reduced
    char *argv0[] = {NULL, argv[4]};
    cartstatehandle = cartstate_new(2, argv0);
  }
  else {
    /*support for non-reduced PDBs*/
    char *argv0[] = {NULL, NULL};
    cartstatehandle = cartstate_new(2, argv0);

    CartState &cs = cartstate_get(cartstatehandle);
    cs.nlig = argc - 4;
    Coor *xx = (Coor *) &cs.x[0];
    int pos = 0;
    cs.ieins[0] = 0;
    for (i = 0; i < cs.nlig; i++) {
      FILE *fil = fopen(argv[i+4], "r");
      Coor *coor;
      read_pdb2(fil,coor,pdbstrings[i],pdblayout[i],cs.natom[i],linecounter[i]);
      if (cs.natom[i]) {
        memcpy(xx,coor,cs.natom[i]*sizeof(Coor));
        delete [] coor;      
        xx += cs.natom[i];	
	pos += cs.natom[i];
      }
      cs.ieins[i] = pos;
    }
    cs.nall = pos;
    memcpy(cs.xori0,cs.x,TOTMAXATOM*3*sizeof(double));
  }

  CartState &cs = cartstate_get(cartstatehandle);  
  memcpy(cs.xori,cs.xori0,TOTMAXATOM*3*sizeof(double));  
  if (modefile != NULL) {
    const int multi = 1;
    read_hm_(modefile,"ligand",cs.nlig, cs.natom, cs.nhm, cs.val, (double *) cs.eig, multi, strlen(modefile), strlen("ligand"));
  }
      
  for (int n = 0; n < enscount; n++) {
    read_ens(cartstatehandle, ens_ligands[n]-1, ens_files[n], 0, morphing[n]);
  }      
      
  //retrieve the parameters needed to read the DOFs
  int nlig; int *nhm;
  cartstate_get_nlig_nhm_(cartstatehandle, nlig,nhm);
  
  int nrdof = 6 * nlig;
  for (int n = 0; n < nlig; n++) nrdof += nhm[n];
  if (nrdof > MAXDOF) {
    fprintf(stderr, "Too many DOFs: %d, MAXDOF=%d\n",nrdof, MAXDOF); 
    exit(1);
  }
  
  //read DOFs and set pivots
  //fpivot contains any pivots read from the DOF file
  double fpivot[3][MAXLIG];
  int auto_pivot, centered_receptor, centered_ligands;  
  int line;
  FILE *fil = read_dof_init_(argv[3], nlig, line, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[3]));
  if (auto_pivot) cartstate_pivot_auto_(cartstatehandle);
  else cartstate_set_pivot_(cartstatehandle, fpivot);  
  
  double *pivot; //the actual pivots (from file or auto-calculated)
  
  int *kai;char4 *tyi;char4 *rgi; int *iei; 
  double *x; int *iaci; double *xlai; int *icop; double *we; int *ieins;
  double *eig; double *xb; double *dmmy1; double *dmmy2;
    
  //get the Cartesian parameters we need for rotation and deformation
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, ieins, eig, pivot, xb, x, dmmy1, dmmy2);
  int *nrens; //the ensemble size for each ligand
  cartstate_get_nrens_(cartstatehandle,nrens);
    
  //get the Cartesian parameters we need for PDB writeout
  cartstate_f_write_pdb_(cartstatehandle,
   nlig, kai,tyi,rgi,iei,x,iaci,xlai,icop,we,ieins); 
 
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
   
  //main loop
  int nstruc = 0;  
  while (1) {
    
    int result = read_dof_(fil, line, nstruc, argv[3], ens, phi, ssi, rot, 
     xa, ya, za, locrests, 
     morph, dlig, nlig, nhm, nrens, morphing, has_locrests,
     seed, label, 0, strlen(argv[3])
    );
    if (result != 0) break;

    if (centered_receptor) { //...then subtract pivot from receptor
      xa[0] -= pivot[0];
      ya[0] -= pivot[MAXLIG];
      za[0] -= pivot[2*MAXLIG];
    }
    

    if (centered_ligands) { //...then subtract pivot from all (other) ligands 
      for (i = 1; i < nlig; i++) {
        xa[i] -= pivot[i];
	ya[i] -= pivot[MAXLIG+i];
	za[i] -= pivot[2*MAXLIG+i];
      }
    }          
    for (i = 0; i < nlig; i++) {

      //Get ensemble differences
      double *ensdp;
      double cmorph;
      double *cmorphdp;
      cartstate_get_ensd_(cartstatehandle, i, ens[i], ensdp,
      morph[i],cmorph, cmorphdp);

      //Apply harmonic modes
      double (&dligp)[MAXMODE] = dlig[i];
      deform_(MAXLIG, 3*MAXATOM, 3*TOTMAXATOM, MAXATOM,MAXMODE, 
        ens[i], ensdp, cmorph, cmorphdp, dligp, nhm, i, ieins, eig, xb, x, dmmy1, dmmy2, 1);
     
      //Compute rotation matrix
      double rotmat[9];
      euler2rotmat_(phi[i],ssi[i],rot[i],rotmat);
      
      //Apply rotation matrix and translation vector
      rotate_(MAXLIG,3*MAXATOM,rotmat,
      xa[i],ya[i],za[i],
      pivot,i,ieins,x);

    }

    CartState &cs = cartstate_get(cartstatehandle);
    Coor *pdb = (Coor *) &cs.x[0];

    int nratoms = cs.nall;
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
     
    /* 
    printf("CONTOUR START\n");
    for (int z = 1; z < extz-1; z++) {
      for (int y = 1; y < g_exty-1; y++) {
        for (int x = 1; x < g_extx-1; x++) {      
          int n = g_exty*g_extx*z + g_extx * y + x;
          print_contour(x,y,z,mapdata_x, mappdb_x, n, threshold);
          print_contour(x,y,z,mapdata_y, mappdb_y, n, threshold);
          print_contour(x,y,z,mapdata_z, mappdb_z, n, threshold);
        }
      }
    }
    printf("CONTOUR END\n");
    */

    //double r = get_corr(vox_map, vox_pdb, nvox_filtered);
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
  
}
