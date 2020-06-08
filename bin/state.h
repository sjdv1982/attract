#ifndef STATE_H
#define STATE_H

#include "max.h"
#include "grid.h"

struct AxSymmetry {
  int ligand;
  int symtype; //0 for ncsym
  double angle; //only for ncsym
  double axis[3];
  double origin[3];
};

struct SymTrans {
  int ligand;
  int targetligand;
  double rotmatsym[9];
  double origin[3];
};

struct CartState {
  /* atom limits */
  int nlig0; //number of ligands before axsym application
  int nlig; //number of ligands after axsym application
  
  int natom[MAXLIG], n3atom[MAXLIG], nres[MAXLIG];
  int ieins[MAXLIG], ieins3[MAXLIG];
  int nall, nall3;

  /* atom values */
  int kai[TOTMAXATOM];      //atom index; not used in calculations
  char4 tyi[TOTMAXATOM];  //atom code; not used in calculations
  char4 rgi[TOTMAXATOM];  //residue code; not used in calculations
  int iei[TOTMAXATOM];      //residue nr
  double x[3*TOTMAXATOM];   //coordinates
  double f[3*TOTMAXATOM];   //forces
  double pivot[3][MAXLIG];  //rotation pivot
  double xb[3*TOTMAXATOM];  //pivot-centered coordinates
  double xori0[3*TOTMAXATOM];  //original coordinates
  double xori[3*TOTMAXATOM];  //original coordinates, after deformation
  int iaci[TOTMAXATOM];     //bead code
  int iaci_old[TOTMAXATOM];     //original bead code; only used in EM data
  double xlai[TOTMAXATOM];  //charge; not used directly in calculations
  int icop[TOTMAXATOM];    //copy counter
  double we[TOTMAXATOM];   //occupancy/weight (current, after lambda)
  double we0[TOTMAXATOM];   //occupancy/weight (original)
  double chai[TOTMAXATOM];  //charge * occupancy
  int use_softcore;		//use softcore potential if set to 1, zero otherwise
  double softcore;	    
  /* modes */
  int nhm[MAXLIG];          //number of modes per ligand
  double val[MAXLIG][MAXMODE];  //force constant per mode
  double eig[MAXLIG][MAXMODE][3*MAXATOM];  //modes

  /* index modes */
  int nihm[MAXLIG];        //number of index modes per ligand
  int index_eig[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE]; //index modes: position of nonzero mode entries
  double index_val[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE]; //index modes: values of nonzero mode entries

  /* copies */
  int ncop[TOTMAXRES][21][11];
  int nmaxco[TOTMAXRES];
  int natco[TOTMAXRES];

  /* forcefield parameters */
  Parameters rbc;  
  Parameters rc;
  Parameters ac;
  Parameters emin;
  Parameters rmin2;
  iParameters ipon;
  iParameters haspar;
  int potshape; //potential shape
  int cdie; //use constant dielectric
  double epsilon; //dielectric constant
  float swi_on; //start (A) of switching
  float swi_off;  //end (A) of switching

  /* grid representations */
  Grid *grids[MAXLIG];  
  
  /* translation table */
  int *transtable;
  
  /*ensemble*/
  int nrens[MAXLIG];
  double *ensd[MAXLIG][MAXENS];   //ensemble delta coordinates
  double **ensw[MAXLIG];   //ensemble pairwise RMSDs (used in MC, only initialized when needed)
  
  /*morphing*/
  int morphing[MAXLIG];
  double *morphd[MAXLIG][MAXENS];   //morphing inter-delta coordinates
  double morph_fconstant;
  
  /*lambda morphing*/
  bool use_lambda;
  double lambda;
  
  /*symmetries*/
  int nsym;
  int symtypes[MAXLIG];
  int sym[MAXLIG][MAXLIG];

  /*axis and non-crystallographic symmetries*/
  int nr_axsyms;
  AxSymmetry axsyms[MAXLIG];
  int nr_symcopies[MAXLIG];
  int symcopies[MAXLIG][24*MAXLIG];
  int nr_symtrans;
  SymTrans symtrans[24*MAXLIG];
  double forcefactor[MAXLIG]; //if 0: no forcerotation
  double forcerotation[MAXLIG][9]; //initialized to the identity matrix
  

  /*location restraints*/
  int has_locrests[MAXLIG];

};

struct MolPair {
  int receptor;
  int ligand;
  int *iactr; //iactr[MAXATOM];
  int *iactl; //iactl[MAXATOM];
  int nonp;              //pairlist: number of pairs
  int *nonr; //nonr[MAXMOLPAIR]  //pairlist: receptor atom
  int *nonl; //nonl[MAXMOLPAIR]  //pairlist: ligand atom
  Grid *grid;
  int pairgen_done;
  int use_energy; //Only calculate energy for the first molpair
};

struct Restraint {
  int type;
  int *selection1; 
  int s1;
  int *selection2; 
  int s2;
  double par1;
  double par2;
  double par3;
  double par4;
  double par5;
  double par6;
  char position_type; //3-bit mask for xyz
  int maxindex;
};

struct MiniState {
  MolPair *pairs;
  int npairs;
  int imc;     //Monte Carlo mode: 0 = off (minfor), 1 = (monte) 2 = (mc_min)
  double mctemp; //Monte Carlo temperature (in KT)
  double mcmtemp; //MCM temperature (in KT)
  double mcscalerot; //Monte Carlo rotation step size (in radians)
  double mcscalecenter; //Monte Carlo translation step size (in A)
  double mcscalemode; //Monte Carlo mode step size (in mode A)
  double mcensprob; //Monte Carlo probability of switching ensemble copies
  int iscore;  //scoring mode: 0 = normal, 1 = scoring, 2 = trajectory
  int ivmax; //max steps
  int imcmax; //max MC steps
  int iori;  //enable orientations
  int itra;  //enable translations
  int ieig;  //enable mode displacement
  int iindex; //enable index mode displacement
  int irst;  //enable CoM restraints
  int fixre; //fix receptor
  double rcut;  //square of distance cutoff
  Restraint *restraints;
  int nr_restraints;  
  int has_globalenergy; //1 = the energy has a global component, globalenergy must be called
  int gravity; //gravity modes: 0 = off, 1 = to global origin, 2 = to receptor origin, 3 = to all other centers;
  double rstk; //gravity force constant
  double restweight; //restraint weight
  bool ghost;
  bool ghost_ligands; //if enabled, ligands don't see each other, only the receptor
};

typedef int (&intarr)[TOTMAXATOM];
typedef double (&dbl3arr)[3*TOTMAXATOM];
typedef double (&dblarr)[TOTMAXATOM];
typedef char (&codearr)[TOTMAXATOM][4];
typedef int (&ncopligtype)[TOTMAXRES][21][11];
typedef int (&copyarr)[TOTMAXRES];
typedef int (&limitarr)[MAXLIG];


extern "C" void read_hm_(const char *hmfile_, const char *hmword_, const int &nlig, const int *natom, int *nhm, double (&vall)[MAXLIG][MAXMODE], double *eigl, const int &multi, int hmfile_len, int hmword_len);

extern "C" void read_indexmode_(const char *hmfile_, const char *hmword_, const int &nlig, int *nhm, int *eigl, double *eigl_val, const int &multi, int hmfile_len, int hmword_len);

extern "C" void read_one_pdb_(
   const int &maxlig, const int &totmaxres, const int &totmaxatom,
   const int &maxatom, const char *pdbfile,
   intarr kai,codearr tyi,codearr rgi,
   intarr iei, dbl3arr x,
   intarr iaci,dblarr xlai,
   intarr icop, dblarr we,dblarr we0,dblarr chai,
   ncopligtype ncoplig, copyarr nmaxco,copyarr natco,
   int &nlig, limitarr nres, limitarr natom, limitarr n3atom,
   int &nall,int &nall3,
   limitarr ieins,limitarr ieins3,
   int pdbfile_len
  );

extern "C" void read_single_pdb_(
   const int &maxlig, const int &totmaxres, const int &totmaxatom,
   const int &maxatom, const char *pdbfile,
   intarr kai,codearr tyi,codearr rgi,
   intarr iei, dbl3arr x,
   intarr iaci,dblarr xlai,
   intarr icop, dblarr we,dblarr we0,dblarr chai,
   ncopligtype ncoplig, copyarr nmaxco,copyarr natco,
   int &nlig, limitarr nres, limitarr natom, limitarr n3atom,
   int &nall,int &nall3,
   limitarr ieins,limitarr ieins3,
   int &i, int &irs,
   int pdbfile_len
  );

extern "C" void read_two_pdbs_(
   const int &maxlig, const int &totmaxres, const int &totmaxatom,
   const int &maxatom, const char *pdbfile1,const char *pdbfile2,
   intarr kai,codearr tyi,codearr rgi,
   intarr iei, dbl3arr x,
   intarr iaci,dblarr xlai,
   intarr icop, dblarr we,dblarr we0,dblarr chai,
   ncopligtype ncoplig, copyarr nmaxco,copyarr natco,
   limitarr nres, limitarr natom, limitarr n3atom,
   int &nall,int &nall3,
   limitarr ieins,limitarr ieins3,
   int pdbfile1_len, int pdbfile2_len
  );

extern "C" void apply_permi_(
  const int &totmaxatom,
  const int &nall,
  dblarr chai,
  const double &permi
);  
#endif
