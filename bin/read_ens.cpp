#include "state.h"
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <cstring>

extern bool exists(const char *f);
extern CartState *cartstates[100];
extern int cartstatesize;

extern CartState &cartstate_get(int handle);

extern "C" void cartstate_select_ligand_(
 const int &handle,
 const int &ligand,
 int &natom,
 int &nres,
 int *&iei,
 double *&x,
 double *&f,
 double *pivot,
 int *&iaci,
 int *&icop,
 double *&we,
 double *&chai,
 int *&ncop,
 int *&nmaxco,
 int *&natco,
 const int &ori
);

extern void read_pdb2(
  FILE *fil, Coor *&x, 
  char **&pdbstrings, bool *&pdblayout,
  int &coorcounter, int &linecounter
);

void read_ens(int cartstatehandle, int ligand, char *ensfile, bool strict, bool morphing) {

  CartState &s = cartstate_get(cartstatehandle);
  CartState *s0 = new CartState;
  CartState &s2 = *s0;
  cartstates[cartstatesize] = s0;
  int cartstatehandle2 = cartstatesize+9990; 

  int natom;
  int nres;
  int *iei;
  double *x;
  double *f;
  double pivot[3];
  int *iaci;
  int *icop;
  double *we;
  double *chai;
  int *ncop;
  int *nmaxco;
  int *natco;
  int ori = 0;

  cartstate_select_ligand_(cartstatehandle, ligand,
   natom,nres,iei,x,f,pivot,iaci,icop,we,chai,ncop,nmaxco,natco,ori);

  char *zeroes;
  if (natom) {
    zeroes = new char[natom*sizeof(int)];
    memset(zeroes, 0, natom*sizeof(int));
  }
  if (strncmp( (char *) icop,zeroes,natom*sizeof(int))) { //copy counters
    fprintf(stderr, "Ensemble file %s has multiple conformation copies\n", ensfile);
    exit(0);
  }
  
  if (!exists(ensfile)) {
    fprintf(stderr, "Ensemble file %s does not exist\n", ensfile);
    exit(0);
  }
  
  FILE *ensf = fopen(ensfile, "r");
  int line = 0;
  int nrens = 0;
  while(!feof(ensf)) {
    line++;
    char fn[1000];
    fscanf(ensf,"%999s\n",fn);
    if (!exists(fn)) {
      fprintf(stderr, "Ensemble file %s, line %d: PDB file %s does not exist\n", ensfile, line, fn);
      exit(0);
    }
    if (strict) {    
      int dmmy=0,dmmy2=0;
      read_single_pdb_(
        MAXLIG, TOTMAXRES, TOTMAXATOM, MAXATOM,
        fn,s2.kai,s2.tyi,s2.rgi,s2.iei,s2.x,s2.iaci,s2.xlai,
        s2.icop,s2.we,s2.we0,s2.chai,s2.ncop,s2.nmaxco,s2.natco,
        s2.nlig,s2.nres,s2.natom,s2.n3atom,s2.nall,s2.nall3,s2.ieins,s2.ieins3,
        dmmy, dmmy2,
        strlen(fn)
      );

    }
    else {    
      char **pdbstrings; bool *pdblayout; int linecounter;
      FILE *fil = fopen(fn, "r");    
      Coor *xx = (Coor *) &s2.x;
      Coor *coor;
      read_pdb2(fil,coor,pdbstrings,pdblayout,s2.natom[0],linecounter);
      if (s2.natom[0]) {
        if (s2.natom[0] >= MAXATOM) {
          fprintf(stderr, "MAXATOM exceeded: %d\n", MAXATOM);
          exit(1);
        }
        memcpy(xx,coor,s2.natom[0]*sizeof(Coor));
        delete [] coor;
      }
    }
    int natom2;
    int nres2;
    int *iei2;
    double *x2;
    double *f2;
    double pivot2[3];
    int *iaci2;
    int *icop2;
    double *we2;
    double *chai2;
    int *ncop2;
    int *nmaxco2;
    int *natco2;
    int ori2 = 0;

    cartstate_select_ligand_(cartstatehandle2, 0,
     natom2,nres2,iei2,x2,f2,pivot2,iaci2,
     icop2,we2,chai2,ncop2,nmaxco2,natco2,ori2);

    if (natom != natom2) { //number of atoms
     if (natom == 0) { //zero atoms, ignore this ligand
       nrens++;
       continue;
     }
     else{
       fprintf(stderr, "Ensemble file %s, line %d: PDB file %s has the wrong number of atoms (%d, expected %d) \n",ensfile, line, fn, natom2, natom);
       exit(0);
     }  
    }
    if (strict) {
      if (nres != nres2) { //number of residues
       fprintf(stderr, "Ensemble file %s, line %d: PDB file %s has the wrong number of residues (%d, expected %d) \n", ensfile, line, fn, nres2, nres);
       exit(0);      
      }
      int d = (iei2[0] - iei[0]);
      for (int n=0; n < natom; n++) iei2[n] -= d;
      if (strncmp( (char *) iei,(char *) iei2,natom*sizeof(int))) { //residue numbers
       fprintf(stderr, "Ensemble file %s, line %d: PDB file %s has a different residue numbering scheme\n", ensfile, line, fn);
       exit(0);      
      }  
      if (strncmp( (char *) iaci,(char *) iaci2,natom*sizeof(int))) { //atom types
       fprintf(stderr, "Ensemble file %s, line %d: PDB file %s has different atom types\n", ensfile, line, fn);
       exit(0);      
      }
      if (strncmp( (char *) icop2,zeroes,natom*sizeof(int))) { //copy counters
       fprintf(stderr, "Ensemble file %s, line %d: PDB file %s has multiple conformation copies\n", ensfile, line, fn);
       exit(0);      
      }
    }
    double *ens = new double[3*MAXATOM];
    for (int n = 0; n < 3*natom; n++) {
      ens[n] = x2[n] - x[n];
    } 
    if (nrens >= MAXENS) {
     fprintf(stderr, "Ensemble file %s, line %d: maximum ensemble size (MAXENS) exceeded\n", ensfile, line);
     exit(0);            
    }
    s.ensd[ligand][nrens] = ens;
    if (nrens > 0 && morphing) {
      double *morphd = new double[3*MAXATOM];
      double *ensd0 = s.ensd[ligand][nrens-1];
      for (int n = 0; n < 3*natom; n++) {
        morphd[n] = ens[n] - ensd0[n];
      }          
      s.morphd[ligand][nrens-1] = morphd;
    }    
    nrens++;
  }
  s.nrens[ligand] = nrens;
  delete s0;
}
