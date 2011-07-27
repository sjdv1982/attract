//Calculated ligand RMSD versus reference structure (identity matrix)

//usage: ./lrmsd structures.dat ligand-unbound.pdb ligand-bound.pdb [ligand2-unbound.pdb ligand2-bound.pdb] [...]
//  structures.dat must have been fitted on the receptor
//  the unbound ligand must have been fitted on the bound ligand

#include "max.h"
#include <cmath>
#include "state.h"
#include <cstdio>

extern int cartstate_new(int argc, char *argv[],bool single=0);

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[3][MAXLIG]);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *(&nhm), int *&ieins, double *&eig, double *&pivot, double *&xb, double *&x,double *&xori, double *&xori0);

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nlm));
extern "C" void cartstate_get_pivot_(const int &handle,double *&pivot);

extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, modes2 &dlig, const int &nlig, const int *nhm, int &seed, char *&label, int f_len);

extern "C" void write_pdb_(
  const int &totmaxatom, const int &maxlig, const int &nlig,
  int *kai, char4 *tyi, char4 *rgi, int *iei, double *x,
  int *iaci, double *xlai, int *icop, double *we, int *ieins);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

extern "C" void rotate_(const int &maxlig,const int &max3atom,double (&rotmat)[9],const double &xa,const double &ya,const double &za,
double *pivot,
int &ijk,int *ieins, double *x);

extern "C" void deform_(const int &maxlig,const int &max3atom, 
const int &totmax3atom, const int &maxatom,const int &maxmode,double (&dligp)[MAXMODE],
int *nhm,int &ijk,int *ieins,double *eig,double *xb,double *x,double *xori,double *xori0); 


extern void read_pdb2(
  FILE *fil, Coor *&x, 
  char **&pdbstrings, bool *&pdblayout,
  int &coorcounter, int &linecounter
);

extern void write_pdb2(
  FILE *fil, const Coor *x, 
  char **pdbstrings, const bool *pdblayout,
  int coorcounter, int linecounter
);

CartState &cartstate_get(int handle);

/* DOFs */
static double phi[MAXLIG];
static double ssi[MAXLIG];
static double rot[MAXLIG];
static double xa[MAXLIG];
static double ya[MAXLIG];
static double za[MAXLIG];
static double dlig[MAXLIG][MAXMODE];
static int seed;
static char *label;

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/lrmsd structures.dat ligand-unbound.pdb ligand-bound.pdb [ligand2-unbound.pdb ligand2-bound.pdb] [...] [...] [...]\n");
  exit(1);
}

extern bool exists(const char *f);

double calc_rmsd(Coor *x, Coor *xrefe, int natoms) {
  double rmsd = 0;
  for (int ii = 0; ii < natoms; ii++) {
    double dx = x[ii][0] - xrefe[ii][0];
    double dy = x[ii][1] - xrefe[ii][1];
    double dz = x[ii][2] - xrefe[ii][2];
    double dsq = dx*dx+dy*dy+dz*dz;
    rmsd += dsq;
  }
  if (rmsd > 0 && natoms > 0) rmsd = sqrt(rmsd/natoms);
  return rmsd;
}

int read_ligands(CartState &cs, char **argv, Coor **bound, Coor *allbound) {
  int ret = 0;
  Coor *xx = (Coor *) &cs.x[0];
  int pos = 0;
  cs.ieins[0] = 0;
  cs.natom[0] = 0;    
  for (int i = 1; i < cs.nlig; i++) {
    char *ub = argv[2*i];
    char *b = argv[2*i+1];
    FILE *fil;
    Coor *coor1; Coor *coor2;
    int ncoor1, ncoor2;
    char **pdbstrings;
    bool *pdblayout;
    int linecounter;

    fil = fopen(ub, "r");
    read_pdb2(fil,coor1,pdbstrings,pdblayout,ncoor1,linecounter);

    if (!ncoor1) {
      fprintf(stderr, "Unbound PDB %s contains no atoms\n", ub);
      exit(1);
    }

    fil = fopen(b, "r");
    read_pdb2(fil,coor2,pdbstrings,pdblayout,ncoor2,linecounter);

    if (!ncoor2) {
      fprintf(stderr, "Bound PDB %s contains no atoms\n", b);
      exit(1);
    }

    if (ncoor1 != ncoor2) {
      fprintf(stderr, "Unbound PDB %s (%d) and bound PDB %s (%d) contains different numbers of atoms\n", ub, ncoor1, b, ncoor2);
      exit(1);
    }

    double rmsd = calc_rmsd(coor1, coor2, ncoor1);
    if (rmsd > 5.0) {
      fprintf(stderr, "Warning: Unbound PDB %s and bound PDB %s have an RMSD of %.3f A, have they been fitted?\n", ub, b, rmsd);    
    }
    ret += ncoor1;
    
    memcpy(xx,coor1,ncoor1*sizeof(Coor));
    memcpy(allbound[pos],coor2,ncoor2*sizeof(Coor));
    delete [] coor1;      
    bound[i] = coor2;

    cs.natom[i] = ncoor1;    
    xx += cs.natom[i];	
    pos += cs.natom[i];        
    cs.ieins[i] = pos;
  }
  memcpy(cs.xori0,cs.x,TOTMAXATOM*3*sizeof(double));
  return ret;
}

int main(int argc, char **argv) {
  int n, i;
  if (argc < 3) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  
  char *modefile = NULL;
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

  if (argc % 2) {
    fprintf(stderr, "Number of arguments must be even: supply unbound ligand, then bound ligand PDB\n");
    exit(1);  
  }
  for (i = 1; i < argc; i++) {
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
  }

  //load the Cartesian parameters and get a handle to it
  int cartstatehandle;
  /*support for non-reduced PDBs*/
  char *argv0[] = {NULL, NULL};
  cartstatehandle = cartstate_new(2, argv0);

  CartState &cs = cartstate_get(cartstatehandle);
  Coor *bound[MAXLIG];
  
  cs.nlig = (argc - 2)/2 + 1;
  Coor allbound[TOTMAXATOM];
  int natoms = read_ligands(cs, argv, bound, allbound);
  
  if (modefile != NULL) {
    CartState &cs = cartstate_get(cartstatehandle);
    const int multi = 1;
    read_hm_(modefile,"ligand",cs.nlig, cs.natom, cs.nhm, cs.val, (double *) cs.eig, multi, strlen(modefile), strlen("ligand"));
  }
      
  //retrieve the parameters needed to read the DOFs
  int nlig, nstruc; int *nhm;
  cartstate_get_nlig_nhm_(cartstatehandle, nlig,nhm);
  
  int nrdof = 6 * nlig;
  for (n = 0; n < nlig; n++) nrdof += nhm[n];
  if (nrdof > MAXDOF) {
    fprintf(stderr, "Too many DOFs: %d, MAXDOF=%d\n",nrdof, MAXDOF); 
    exit(1);
  }
  
  //read DOFs and set pivots
  //fpivot contains any pivots read from the DOF file

  double fpivot[3][MAXLIG]; 
  int auto_pivot, centered_receptor, centered_ligands;  
  //read_dof_(argv[1], phi, ssi, rot, xa, ya, za, dlig, nlig, nhm, seed, label, nstruc, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  int line;
  FILE *fil = read_dof_init_(argv[1], nlig, line, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  
  if (auto_pivot) cartstate_pivot_auto_(cartstatehandle);
  else cartstate_set_pivot_(cartstatehandle, fpivot);  
  
  double *pivot; //the actual pivots (from file or auto-calculated)
  cartstate_get_pivot_(cartstatehandle,pivot);

  //get the Cartesian parameters we need for rotation and deformation
  double *x; int *ieins;double *eig; double *xb; double *xori; double *xori0;
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, ieins, eig, pivot, xb, x, xori, xori0);

  nstruc = 0;
  while (1) {
    int result = read_dof_(fil, line, nstruc, argv[1], phi, ssi, rot, xa, ya, za, dlig, nlig, nhm, seed, label, strlen(argv[1]));
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
    if (fabs(phi[0]) > 0.001 ||
        fabs(ssi[0]) > 0.001 ||
	fabs(rot[0]) > 0.001 ||
	fabs(xa[0]) > 0.001 ||
	fabs(ya[0]) > 0.001 ||
	fabs(za[0]) > 0.001) 
    {
      fprintf(stderr, "ERROR: Structures have not yet been fitted\n");
      return 1;
    }
    for (i = 1; i < nlig; i++) {
     
      //Apply harmonic modes
      double (&dligp)[20] = dlig[i];
      deform_(MAXLIG, 3*MAXATOM, 3*TOTMAXATOM,MAXATOM, MAXMODE, 
        dligp, nhm, i, ieins, eig, xb, x,xori, xori0);
     
      //Compute rotation matrix
      double rotmat[9];
      euler2rotmat_(phi[i],ssi[i],rot[i],rotmat);
      
      //Apply rotation matrix and translation vector
      rotate_(MAXLIG,3*MAXATOM,rotmat,
      xa[i],ya[i],za[i],
      pivot,i,ieins,x);
      
    }
    Coor *x2 = (Coor *) x;

    double rmsd = calc_rmsd(x2, allbound, natoms);
    printf("l-RMSD %.3f", rmsd);    
    if (nlig > 2) {
      for (i = 1; i < nlig; i++) {
        rmsd = calc_rmsd(x2+ieins[i-1], allbound+ieins[i-1], ieins[i]-ieins[i-1]);
        printf(" %.3f", rmsd);
      }
    }
    printf("\n");
  }
}
