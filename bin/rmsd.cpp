//Calculated RMSD versus reference structure (identity matrix)

//usage: ./rmsd structures.dat ligand-unbound.pdb ligand-bound.pdb [ligand2-unbound.pdb ligand2-bound.pdb] [...] [--modes <modefile>] [--ens <ligand nr> <ensemble file>] [--locrest <ligand>]
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
  int *(&nhm), int *(&nihm), int *&ieins, double *&eig, int *&index_eig, double *&index_val, double *&pivot, double *&xb, double *&x,double *&xori, double *&xori0);

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nlm), int *(&nihm));
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

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, coors2 &locrests, dof2 &morph, modes2 &dlig, const int &nlig, const int *nhm, const int *nihm, const int *nrens0, const int *morphing, const int *has_locrests, int &seed, char *&label, const int &all_labels, int f_len);

extern "C" void write_pdb_(
  const int &totmaxatom, const int &maxlig, const int &nlig,
  int *kai, char4 *tyi, char4 *rgi, int *iei, double *x,
  int *iaci, double *xlai, int *icop, double *we, int *ieins);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

extern "C" void rotate_(const int &maxlig,const int &max3atom,double (&rotmat)[9],const double &xa,const double &ya,const double &za,
double *pivot,
int &ijk,int *ieins, double *x);

extern "C" void deform_(int &ens, double *ensdp, const double &cmorph, const double *cmorphdp,
double (&dligp)[MAXMODE+MAXINDEXMODE],
int *nhm,int *nihm,int &ijk,int *ieins,double *eig,int *index_eig, double *index_val,double *xb,double *x,double *xori,double *xori0, const int &do_morph);


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
static double dlig[MAXLIG][MAXMODE+MAXINDEXMODE];
static int seed;
static char *label;

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/lrmsd structures.dat ligand-unbound.pdb ligand-bound.pdb [ligand2-unbound.pdb ligand2-bound.pdb] [...] [...] [--modes <modefile>] [--imodes <indexmodefile>] [--ens <ligand nr> <ensemble file>]\n");
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
  for (int i = 0; i < cs.nlig; i++) {
    char *ub = argv[2*i+2];
    char *b = argv[2*i+3];
    printf("%s %s %d\n", b, ub, cs.nlig);
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
    memcpy(&allbound[pos],coor2,ncoor2*sizeof(Coor));
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

int morphing[MAXLIG];
int main(int argc, char **argv) {
  memset(morphing,0,MAXLIG*sizeof(int));
  int n, i;

  coors2 locrests;
  int has_locrests[MAXLIG];
  memset(has_locrests, 0, MAXLIG*sizeof(int));

  if (argc < 3) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  
  char *modefile = NULL;
  char *indexmodefile = NULL;
  int enscount = 0;
  int ens_ligands[MAXLIG];
  char *ens_files[MAXLIG];
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
  for (int n = 1; n < argc-1; n++) {
    if (!strcmp(argv[n],"--imodes")) {
      indexmodefile = argv[n+1];
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+2 < argc) memcpy(argv2+n,argv+n+2,(argc-n-2)*sizeof(char*));
      argv = argv2;
      argc -= 2;
      break;
    }
  }
  for (int n = 1; n < argc-2; n++) {
    if (!strcmp(argv[n],"--ens")) {
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
  
  cs.nlig = (argc - 2)/2;
  Coor allbound[TOTMAXATOM];
  int natoms = read_ligands(cs, argv, bound, allbound);
  
  if (modefile != NULL) {
    CartState &cs = cartstate_get(cartstatehandle);
    const int multi = 1;
    read_hm_(modefile,"ligand",cs.nlig, cs.natom, cs.nhm, cs.val, (double *) cs.eig, multi, strlen(modefile), strlen("ligand"));
  }
  if (indexmodefile != NULL) {
	    CartState &cs = cartstate_get(cartstatehandle);
	    const int multi = 1;
	    read_indexmode_(indexmodefile,"ligand",cs.nlig, cs.nihm, (int *) cs.index_eig, (double *) cs.index_val, multi, strlen(indexmodefile), strlen("ligand"));
  }
  for (int n = 0; n < enscount; n++) {
    read_ens(cartstatehandle, ens_ligands[n]-1, ens_files[n], 0, 0);
  }      
      
  //retrieve the parameters needed to read the DOFs
  int nlig, nstruc; int *nhm; int *nihm;
  cartstate_get_nlig_nhm_(cartstatehandle, nlig,nhm, nihm);
  
  int nrdof = 6 * nlig;
  for (n = 0; n < nlig; n++) {
	  nrdof += nhm[n];
	  nrdof += nihm[n];
  }
  if (nrdof > MAXDOF) {
    fprintf(stderr, "Too many DOFs: %d, MAXDOF=%d\n",nrdof, MAXDOF); 
    exit(1);
  }
  
  //read DOFs and set pivots
  //fpivot contains any pivots read from the DOF file

  double fpivot[3][MAXLIG]; 
  int auto_pivot, centered_receptor, centered_ligands;  
  int line;
  FILE *fil = read_dof_init_(argv[1], nlig, line, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  
  if (auto_pivot) cartstate_pivot_auto_(cartstatehandle);
  else cartstate_set_pivot_(cartstatehandle, fpivot);  
  
  double *pivot; //the actual pivots (from file or auto-calculated)
  cartstate_get_pivot_(cartstatehandle,pivot);
 
  //get the Cartesian parameters we need for rotation and deformation
  double *x; int *ieins;double *eig; int *index_eig; double *index_val; double *xb; double *xori; double *xori0;
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, nihm, ieins, eig, index_eig, index_val,pivot, xb, x, xori, xori0);
  int *nrens; //the ensemble size for each ligand
  cartstate_get_nrens_(cartstatehandle,nrens);

  nstruc = 0;

  while (1) {
    int result = read_dof_(fil, line, nstruc, argv[1], ens, phi, ssi, rot, 
     xa, ya, za, locrests, 
     morph, dlig, nlig, nhm, nihm, nrens, morphing, has_locrests,
     seed, label, 0, strlen(argv[1])
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
      double dmmy1; double *dmmy2;      
      cartstate_get_ensd_(cartstatehandle, i, ens[i], ensdp, -1, dmmy1, dmmy2);
     
      //Apply harmonic modes and ensemble differences
      double (&dligp)[MAXMODE+MAXINDEXMODE] = dlig[i];
      deform_(ens[i], ensdp, -1, NULL, dligp, nhm, nihm, i, ieins, eig, index_eig, index_val, xb, x,xori,xori0,0);
     
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
    printf("RMSD %.3f", rmsd);    
    if (nlig > 1) {
      for (i = 0; i < nlig; i++) {
        rmsd = calc_rmsd(x2+ieins[i-1], allbound+ieins[i-1], ieins[i]-ieins[i-1]);
        printf(" %.3f", rmsd);
      }
    }
    printf("\n");
  }
}
