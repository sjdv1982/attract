//Calculated ligand RMSD versus reference structure (identity matrix)
//fix the receptor first, subtract centers, and do any other operation
// so that the identity matrix indeed is the reference

//usage: ./lrmsd structures.dat receptor.pdb [ligand.pdb]
//  if no ligand.pdb, receptor.pdb is a multi-ligand PDB file 


#include "max.h"
#include <cmath>

extern int cartstate_new(int argc, char *argv[],bool single=0);

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[3][MAXLIG]);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_f_write_pdb_(
  const int &handle,
  int &nlig, int *&kai, char4 *&tyi, char4 *&rgi, int *&iei, double *&x,
  int *&iaci, double *&xlai, int *&icop, double *&we, int *&ieins);

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *(&nhm), int *&ieins, double *&eig, double *&pivot, double *&xb, double *&x,double *&xori, double *&xori0);

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nlm));
extern "C" void cartstate_get_pivot_(const int &handle,double *&pivot);

extern "C" void read_dofs_(const char *f_, dof &phi, dof &ssi, dof &rot, dof &xa, dof &ya, dof &za, modes &dlig, const int &nlig, const int *nhm, int (&seed)[MAXSTRUC], char *label[MAXSTRUC], int &nstruc, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

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

/*support for non-reduced PDBs*/
#include "state.h"
#include <cstdio>

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
static double phi[MAXSTRUC][MAXLIG];
static double ssi[MAXSTRUC][MAXLIG];
static double rot[MAXSTRUC][MAXLIG];
static double xa[MAXSTRUC][MAXLIG];
static double ya[MAXSTRUC][MAXLIG];
static double za[MAXSTRUC][MAXLIG];
static double dlig[MAXSTRUC][MAXLIG][MAXMODE];

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/lrmsd structures.dat receptor.pdb [ligand.pdb] [...] [...]\n");
  exit(1);
}

extern bool exists(const char *f);

int main(int argc, char *argv[]) {
  int n, i;
  if (argc < 3) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  for (i = 1; i <= 2; i++) {
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
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
  char **pdbstrings[MAXLIG]; bool *pdblayout[MAXLIG]; int linecounter[MAXLIG];

  //load the Cartesian parameters and get a handle to it
  int cartstatehandle;
  if (argc == 3) { //one PDB, reduced or non-reduced
    char *argv0[] = {NULL, argv[2]};
    cartstatehandle = cartstate_new(2, argv0);
  }
  else {
    //char *argv0[] = {NULL, argv[2], argv[3]};
    //cartstatehandle = cartstate_new(3, argv0);

    /*support for non-reduced PDBs*/
    char *argv0[] = {NULL, NULL};
    cartstatehandle = cartstate_new(2, argv0);

    CartState &cs = cartstate_get(cartstatehandle);
    cs.nlig = argc - 2;
    Coor *xx = (Coor *) &cs.x[0];
    int pos = 0;
    cs.ieins[0] = 0;
    for (i = 0; i < cs.nlig; i++) {
      FILE *fil = fopen(argv[i+2], "r");
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
    memcpy(cs.xori0,cs.x,TOTMAXATOM*3*sizeof(double));
  }
  
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
  int seed[MAXSTRUC]; char *label[MAXSTRUC];
  int auto_pivot, centered_receptor, centered_ligands;  
  read_dofs_(argv[1], phi, ssi, rot, xa, ya, za, dlig, nlig, nhm, seed, label, nstruc, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  if (auto_pivot) cartstate_pivot_auto_(cartstatehandle);
  else cartstate_set_pivot_(cartstatehandle, fpivot);  
  
  double *pivot; //the actual pivots (from file or auto-calculated)
  
  int *kai;char4 *tyi;char4 *rgi; int *iei; 
  double *x; int *iaci; double *xlai; int *icop; double *we; int *ieins;
  double *eig; double *xb; double *dmmy1; double *dmmy2;
    
  //get the Cartesian parameters we need for rotation and deformation
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, ieins, eig, pivot, xb, x, dmmy1, dmmy2);
    
  if (centered_receptor) { //...then subtract pivot from receptor
    cartstate_get_pivot_(cartstatehandle,pivot);
    for (n = 0; n < nstruc; n++) {
      xa[n][0] -= pivot[0];
      ya[n][0] -= pivot[MAXLIG];
      za[n][0] -= pivot[2*MAXLIG];
    }
  }  

  double xaref=0,yaref=0,zaref=0;
  if (centered_ligands) { //...then subtract pivot from all (other) ligands 
    cartstate_get_pivot_(cartstatehandle,pivot);
    for (n = 0; n < nstruc; n++) {
      for (i = 1; i < nlig; i++) {
        xa[n][i] -= pivot[i];
	ya[n][i] -= pivot[MAXLIG+i];
	za[n][i] -= pivot[2*MAXLIG+i];
      }
      xaref -= pivot[i];
      yaref -= pivot[MAXLIG+i];
      zaref -= pivot[2*MAXLIG+i];
    }
  }  
          

  //get the Cartesian parameters we need for PDB writeout
  cartstate_f_write_pdb_(cartstatehandle,
   nlig, kai,tyi,rgi,iei,x,iaci,xlai,icop,we,ieins); 
  
  //main loop  

  for (n = 0; n < nstruc; n++) {    
    if (fabs(phi[n][0]) > 0.001 ||
        fabs(ssi[n][0]) > 0.001 ||
	fabs(rot[n][0]) > 0.001 ||
	fabs(xa[n][0]) > 0.001 ||
	fabs(ya[n][0]) > 0.001 ||
	fabs(za[n][0]) > 0.001) 
    {
      fprintf(stderr, "ERROR: Structures have not yet been fitted\n");
      return 1;
    }
  }
 
  //Calculate reference
  for (i = 1; i < nlig; i++) {
      
    //Apply harmonic modes
    double (&dligp)[20] = dlig[n][i];
    deform_(MAXLIG, 3*MAXATOM, 3*TOTMAXATOM,MAXATOM, MAXMODE, 
      dligp, nhm, i, ieins, eig, xb, x,dmmy1,dmmy2);

    //Compute rotation matrix
    double rotmat[9];
    euler2rotmat_(0,0,0,rotmat);

    //Apply rotation matrix and translation vector
    rotate_(MAXLIG,3*MAXATOM,rotmat,
    xaref,yaref,zaref,
    pivot,i,ieins,x);
  }
  int natoms = ieins[nlig-1];
  double *xrefe = new double[3*natoms];
  memcpy(xrefe,x,3*natoms*sizeof(double));

  for (n = 0; n < nstruc; n++) {    
    for (i = 1; i < nlig; i++) {
      
      //Apply harmonic modes
      double (&dligp)[20] = dlig[n][i];
      deform_(MAXLIG, 3*MAXATOM, 3*TOTMAXATOM,MAXATOM, MAXMODE, 
        dligp, nhm, i, ieins, eig, xb, x,dmmy1,dmmy2);
     
      //Compute rotation matrix
      double rotmat[9];
      euler2rotmat_(phi[n][i],ssi[n][i],rot[n][i],rotmat);
      
      //Apply rotation matrix and translation vector
      rotate_(MAXLIG,3*MAXATOM,rotmat,
      xa[n][i],ya[n][i],za[n][i],
      pivot,i,ieins,x);
      
    }
    double rmsd = 0;
    for (int ii = ieins[0]; ii < natoms; ii++) {
      double dx = x[3*ii] - xrefe[3*ii];
      double dy = x[3*ii+1] - xrefe[3*ii+1];
      double dz = x[3*ii+2] - xrefe[3*ii+2];
      double dsq = dx*dx+dy*dy+dz*dz;
      rmsd += dsq;
    }
    int cnatoms = natoms - ieins[0];
    if (rmsd > 0 && cnatoms > 0) rmsd = sqrt(rmsd/cnatoms);
    printf("l-RMSD %.3f", rmsd);
    if (nlig > 2) {
      for (i = 1; i < nlig; i++) {
	rmsd = 0;
	for (int ii = ieins[i-1]; ii < ieins[i]; ii++) {
	  double dx = x[3*ii] - xrefe[3*ii];
	  double dy = x[3*ii+1] - xrefe[3*ii+1];
	  double dz = x[3*ii+2] - xrefe[3*ii+2];
	  double dsq = dx*dx+dy*dy+dz*dz;
	  rmsd += dsq;
	}
	int cnatoms = ieins[i]-ieins[i-1];
	if (rmsd > 0 && cnatoms > 0) rmsd = sqrt(rmsd/cnatoms);
	printf(" %.3f", rmsd);
      }
    }
    printf("\n");    
  }  
}
