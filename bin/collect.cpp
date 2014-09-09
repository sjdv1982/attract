//Converts DOF to PDB
//does not eliminate redundant solutions

//usage: ./collect structures.dat receptor.pdb [ligand.pdb] [...] [...] [--modes <modefile>] [--ens/--morph <ligand nr> <ensemble file>]
//  if no ligand.pdb, receptor.pdb is a multi-ligand PDB file 


#include "max.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream> // for ostringstream
#include <string>

extern int cartstate_new(int argc, char *argv[],bool single=0);

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[3][MAXLIG]);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_f_write_pdb_(
  const int &handle,
  int &nlig, int *&kai, char4 *&tyi, char4 *&rgi, int *&iei, double *&x,
  int *&iaci, double *&xlai, int *&icop, double *&we, int *&ieins);

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *(&nhm), int *(&nihm), int *&ieins, double *&eig, int *&index_eig, double *&index_val,
  double *&pivot, double *&xb, double *&x,double *&xori, double *&xori0);

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
int *nhm, int *nihm,int &ijk,int *ieins,double *eig,int *index_eig, double *index_val,double *xb,double *x,double *xori,double *xori0, const int &do_morph);

/*support for non-reduced PDBs*/
#include "state.h"

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
  fprintf(stderr, "usage: collect structures.dat receptor.pdb [...] [...] [--modes <modefile>] [--imodes <indexmodefile>] [--ens/--morph <ligand nr> <ensemble file>] [--name <iATTRACT naming scheme>]");
  exit(1);
}

extern bool exists(const char *f);

int main(int argc, char *argv[]) {
  int i;
  if (argc < 3) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  for (i = 1; i <= 2; i++) {
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
  }

  coors2 locrests;
  int has_locrests[MAXLIG];
  memset(has_locrests, 0, MAXLIG*sizeof(int));
  
  char *modefile = NULL;
  char *indexmodefile = NULL;
  char *name = NULL;
  int enscount = 0;
  int ens_ligands[MAXLIG];
  char *ens_files[MAXLIG];
  int morphing[MAXLIG];
  memset(morphing,0,MAXLIG*sizeof(int));
  bool single = 0;

  for (int n = 1; n < argc; n++) {
    if (!strcmp(argv[n],"--single")) {
      single = 1;
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+1 < argc) memcpy(argv2+n,argv+n+1,(argc-n-1)*sizeof(char*));
      argv = argv2;
      argc -= 1;
      break;
    }
  }
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
  for (int n = 1; n < argc-1; n++) {
    if (!strcmp(argv[n],"--name")) {
      name = argv[n+1];
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


  char **pdbstrings[MAXLIG]; bool *pdblayout[MAXLIG]; int linecounter[MAXLIG];

  //load the Cartesian parameters and get a handle to it
  int cartstatehandle;
  if (argc == 3 && single == 0) { //one PDB, reduced
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
      if (!fil) {
        fprintf(stderr, "File %s does not exist\n", argv[i+2]);
        exit(1);
      }
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

  CartState &cs = cartstate_get(cartstatehandle);  
  memcpy(cs.xori,cs.xori0,TOTMAXATOM*3*sizeof(double));  
  if (modefile != NULL) {
    const int multi = 1;
    read_hm_(modefile,"ligand",cs.nlig, cs.natom, cs.nhm, cs.val, (double *) cs.eig, multi, strlen(modefile), strlen("ligand"));
  }

  if (indexmodefile != NULL) {
	  const int multi = 1;
	      read_indexmode_(indexmodefile,"ligand",cs.nlig, cs.nihm, (int *) cs.index_eig, (double *) cs.index_val, multi, strlen(indexmodefile), strlen("ligand"));
  }
  for (int n = 0; n < enscount; n++) {
    read_ens(cartstatehandle, ens_ligands[n]-1, ens_files[n], 0, morphing[n]);
  }      
  //retrieve the parameters needed to read the DOFs
  int nlig; int *nhm; int *nihm;
  cartstate_get_nlig_nhm_(cartstatehandle, nlig,nhm,nihm);

  int nrdof = 6 * nlig;
  for (int n = 0; n < nlig; n++) {
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
  
  int *kai;char4 *tyi;char4 *rgi; int *iei; 
  double *x; int *iaci; double *xlai; int *icop; double *we; int *ieins;
  double *eig; int *index_eig; double *index_val; double *xb; double *dmmy1; double *dmmy2;

  //get the Cartesian parameters we need for rotation and deformation
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, nihm, ieins, eig, index_eig, index_val, pivot, xb, x, dmmy1, dmmy2);
  int *nrens; //the ensemble size for each ligand
  cartstate_get_nrens_(cartstatehandle,nrens);

  //get the Cartesian parameters we need for PDB writeout
  cartstate_f_write_pdb_(cartstatehandle,
   nlig, kai,tyi,rgi,iei,x,iaci,xlai,icop,we,ieins); 
  
  //main loop
  int nstruc = 0;  
  while (1) {
    if (name != NULL){
       std::ostringstream out; 
       out << "flexm-" << nstruc+1 << name << ".dat";
       std::string save;
       save = out.str();
       const char * newindexmodefile;
       newindexmodefile = save.c_str();
       if (newindexmodefile != NULL) {
	  const int multi = 1;
	  read_indexmode_(newindexmodefile,"ligand",cs.nlig, cs.nihm, (int *) cs.index_eig, (double *) cs.index_val, multi, strlen(newindexmodefile), strlen("ligand"));
       }
    }
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
      double cmorph;
      double *cmorphdp;
      cartstate_get_ensd_(cartstatehandle, i, ens[i], ensdp,
      morph[i],cmorph, cmorphdp);

      //Apply harmonic modes
      double (&dligp)[MAXMODE+MAXINDEXMODE] = dlig[i];
      deform_(ens[i], ensdp, cmorph, cmorphdp, dligp, nhm, nihm, i, ieins, eig,index_eig, index_val,
    		  xb, x, dmmy1, dmmy2, 1);
     
      //Compute rotation matrix
      double rotmat[9];
      euler2rotmat_(phi[i],ssi[i],rot[i],rotmat);
      
      //Apply rotation matrix and translation vector
      rotate_(MAXLIG,3*MAXATOM,rotmat,
      xa[i],ya[i],za[i],
      pivot,i,ieins,x);

    }
    //Write out PDB
    if (argc == 3 && single == 0) {
      write_pdb_(TOTMAXATOM,MAXLIG, nlig,
       kai,tyi,rgi,iei,x,iaci,xlai,icop,we, ieins);
    }
    else {
      printf("MODEL %d\n", nstruc);
      /*support for non-reduced PDBs*/
      CartState &cs = cartstate_get(cartstatehandle);
      Coor *xx = (Coor *) &cs.x[0];
      for (int i = 0; i < cs.nlig; i++) {	  
        write_pdb2(stdout,xx,pdbstrings[i],
         pdblayout[i],cs.natom[i],linecounter[i]);
	xx += cs.natom[i];
	printf("TER\n");      
      }    
      printf("ENDMDL\n");      
    }
  }  
}
