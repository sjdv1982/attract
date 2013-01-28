/*
c main program
c This program performs docking a docking search on a list.
c usage: $path/attract structures.dat parameterfile receptor.pdb [ligand.pdb]
c Author: Martin Zacharias, Jacobs University Bremen
c Adapted by Sjoerd de Vries
*/

#include "max.h"
#include <cmath>
#include <sys/mman.h>
#include <cstdio>

extern "C" void print_struc_(
 const int &seed,
 char *label0, 
 const double &energy,
 const double *energies,
 const int &nlig,
 const int *ens, 
 const double *phi,
 const double *ssi,
 const double *rot,
 const double *xa,
 const double *ya,
 const double *za,
 const coors2 &locrests, 
 const double *morph,
 const int *nhm,
 const modes2 &dlig, 
 const int *has_locrests,
 int len_label
); 

extern int cartstate_new(int argc, char *argv[],bool single=0);
extern "C" void cartstate_translate_atomtypes_(const int &handle);
extern bool exists(const char *);

extern "C" void cartstate_translate_atomtypes_(const int &handle);
extern "C" void cartstate_prepare_axsym_(const int &cartstatehandle);

extern void parse_options(int ministatehandle, int cartstatehandle, int nlig, int argc, char *argv[]);

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[3][MAXLIG]);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_f_write_pdb_(
  const int &handle,
  int &nlig, int *&kai, char4 *&tyi, char4 *&rgi, int *&iei, double *&x,
  int *&iaci, double *&xlai, int *&icop, double *&we, int *&ieins);

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *(&nhm), int *&ieins, double *&eig, double *&pivot, double *&xb, double *&x,double *&xori, double *&xori0);

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nhm));
extern "C" void cartstate_get_pivot_(const int &handle,double *&pivot);
extern "C" void cartstate_get_nrens_(const int &handle,int *&nrens);
extern "C" void cartstate_get_has_locrests_(const int &handle,int *&has_locrests);
extern "C" void cartstate_get_morphing_(const int &handle,int *&morphing);

extern "C" void cartstate_get_ensd_(const int &handle,
  const int &ligand,
  const int &ens,
  double *&ensd,
  const double &morph,
  double &cmorph,
  double *&cmorphd
  );


extern "C" int ministate_new_();

extern "C" void ministate_iscore_imc_(const int &handle, int &iscore, int &imc);

extern "C" void cartstate_apply_epsilon_(const int  &cartstatehandle);

extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, coors2 &locrests, dof2 &morph, modes2 &dlig, const int &nlig, const int *nhm, const int *nrens0, const int *morphing, const int *has_locrests, int &seed, char *&label, const int &all_labels, int f_len);

extern "C" void minfor_(
const int &maxatom, const int &totmaxatom, const int &maxres,const int &totmaxres,
const int &maxlig, const int &maxdof, const int &maxmode,const int &maxmolpair,
const int &cartstatehandle,const int &ministatehandle, 
int *nhm, const int &nlig, 
int *ens, double *phi, double *ssi, double *rot, double *xa, double *ya, double *za, double *morph, 
double *dlig, 
double *locrests, int *has_locrests,
const int &seed, char *label, double &energy, double *energies, int &lablen);

extern "C" void monte_(
const int &maxatom, const int &totmaxatom, const int &maxres,const int &totmaxres,
const int &maxlig, const int &maxdof, const int &maxmode,const int &maxmolpair,
const int &cartstatehandle,const int &ministatehandle, 
int *nhm, const int &nlig, 
int *ens,  double *phi, double *ssi, double *rot, double *xa, double *ya, double *za, double *morph,
double *dlig, 
double *locrests, int *has_locrests,
const int &seed, char *label, double &energy, double *energies, int &lablen);

extern "C" void write_pdb_(
  const int &totmaxatom, const int &maxlig, const int &nlig,
  int *kai, char4 *tyi, char4 *rgi, int *iei, double *x,
  int *iaci, double *xlai, int *icop, double *we, int *ieins);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

extern "C" void rotate_(const int &maxlig,const int &totmax3atom,double (&rotmat)[9],double &xa,double &ya,double &za,
double *pivot,
int &ijk,int *ieins, double *x);

extern "C" void deform_(const int &maxlig,const int &max3atom, 
const int &totmax3atom, const int &maxatom,const int &maxmode,
int &ens, double *ensdp, const double &cmorph, const double *cmorphdp, 
double (&dligp)[MAXMODE], 
int *nhm,int &ijk,int *ieins,double *eig,double *xb,double *x,double *xori,double *xori0, const int &do_morph); 

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
static coors2 locrests;
static int seed;
static char *label;

#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/attract structures.dat parameterfile receptor.pdb [ligand.pdb] [options]\n");
  exit(1);
}

extern "C" void cartstate_get_nlig_nhm_(const int &handle, int &nlig, int *(&nlm));

int main(int argc, char *argv[]) {
  memset(dlig, 0, sizeof(double) * MAXLIG * MAXMODE);
  int n, i;
  int argc0;
  for (argc0 = 1; argc0 < argc; argc0++) {
    if (argv[argc0][0] == '-') break;
  }
  if (argc0 < 4) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  if (argc0 > 5) {
    fprintf(stderr, "Too many arguments\n"); usage();
  }
  for (i = 1; i <= 3; i++) {
    if (!exists(argv[i])) {
      fprintf(stderr, "File %s does not exist\n", argv[i]);
      exit(1);
    }
  }

  //load the Cartesian parameters and get a handle to it
  int cartstatehandle = cartstate_new(argc0-2, argv+2);
      
  //check number of DOFs    
  int nlig; int *nhm;
  cartstate_get_nlig_nhm_(cartstatehandle, nlig,nhm);  
  int nrdof = 6 * nlig;
  for (n = 0; n < nlig; n++) nrdof += nhm[n];
  if (nrdof > MAXDOF) {
    fprintf(stderr, "Too many DOFs: %d, MAXDOF=%d\n",nrdof, MAXDOF); 
    exit(1);
  }
  
  //parse options
  int ministatehandle = ministate_new_();  
  parse_options(ministatehandle, cartstatehandle, nlig, argc-argc0,argv+argc0);
    
  cartstate_translate_atomtypes_(cartstatehandle);
  cartstate_apply_epsilon_(cartstatehandle);
    
  //read DOFs and set pivots
  //fpivot contains any pivots read from the DOF file
  double fpivot[3][MAXLIG]; 
  int auto_pivot, centered_receptor, centered_ligands;  
  int line;
  FILE *fil = read_dof_init_(argv[1], nlig, line, fpivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  
  if (auto_pivot) cartstate_pivot_auto_(cartstatehandle);
  else cartstate_set_pivot_(cartstatehandle, fpivot);  
  
  double *pivot; //the actual pivots (from file or auto-calculated)
  int *nrens; //the ensemble size for each ligand
  int *morphing; //flag: is morphing active or not for each ligand
  int *has_locrests; //flag: are their location restraints or not for each ligand
  cartstate_get_pivot_(cartstatehandle,pivot);
  cartstate_get_nrens_(cartstatehandle,nrens);
  cartstate_get_morphing_(cartstatehandle,morphing);
  cartstate_get_has_locrests_(cartstatehandle,has_locrests);
  
  //get the Cartesian parameters we need for rotation and deformation
  double *x; int *ieins;double *eig; double *xb; double *xori; double *xori0;
  cartstate_f_rotdeform_(cartstatehandle,
   nhm, ieins, eig, pivot, xb, x, xori, xori0);

  int nstruc = 0;
  int iscore, imc;
  ministate_iscore_imc_(ministatehandle, iscore, imc);
  if (iscore != 1) {
    printf("## Command line arguments:");
    for (i = 1; i < argc; i++) {
      printf(" %s", argv[i]);
    }
    printf("\n");
    for (i = 0; i < nlig; i++) {
      printf("#pivot %d %.3f %.3f %.3f\n", 
        i+1, pivot[i], pivot[MAXLIG+i], pivot[2*MAXLIG+i]);
    }              
    printf("#centered receptor: false\n");
    printf("#centered ligands: false\n");
  }
  cartstate_prepare_axsym_(cartstatehandle);
  while (1) {
    int result = read_dof_(fil, line, nstruc, argv[1], ens, phi, ssi, rot, 
     xa, ya, za, locrests, 
     morph, dlig, nlig, nhm, nrens, morphing, has_locrests,
     seed, label, 0, strlen(argv[1])
    );
    if (result != 0) {
      break;
    }  

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

    double energy; double energies[8];
    for (int l = 0; l < nlig; l++) {      

      //Get ensemble differences
      double *ensdp;
      double cmorph;
      double *cmorphdp;
      cartstate_get_ensd_(cartstatehandle, l, ens[l], ensdp,
      morph[l],cmorph, cmorphdp);
                        
      //Apply harmonic modes
      double (&dligp)[MAXMODE] = dlig[l];
      deform_(MAXLIG, 3*MAXATOM, 3*TOTMAXATOM, MAXATOM,MAXMODE, 
        ens[l], ensdp, cmorph, cmorphdp, dligp, nhm, l, ieins, eig, xb, x, xori, xori0, 1);

      //Compute rotation matrix
      double rotmat[9];
      euler2rotmat_(phi[l],ssi[l],rot[l],rotmat);

      //Apply rotation matrix and translation vector
      rotate_(MAXLIG,3*TOTMAXATOM,rotmat,
      xa[l],ya[l],za[l],
      pivot,l,ieins,x);
    }
    int lablen = 1;
    if (label != NULL) lablen = strlen(label);
    if (imc == 0) {
      minfor_(
        MAXATOM,TOTMAXATOM,MAXRES,TOTMAXRES,
        MAXLIG,MAXDOF,MAXMODE,MAXMOLPAIR,
        cartstatehandle, ministatehandle,
        nhm, nlig,
        &ens[0], &phi[0], &ssi[0], &rot[0], 
        &xa[0], &ya[0], &za[0], 
        &morph[0], &dlig[0][0],
        &locrests[0][0], has_locrests,
        seed, label,
        energy, energies, lablen
      );
    }
    else {
      monte_(
        MAXATOM,TOTMAXATOM,MAXRES,TOTMAXRES,
        MAXLIG,MAXDOF,MAXMODE,MAXMOLPAIR,
        cartstatehandle, ministatehandle,
        nhm, nlig,
        &ens[0], &phi[0], &ssi[0], &rot[0], 
        &xa[0], &ya[0], &za[0], 
        &morph[0], &dlig[0][0],
        &locrests[0][0], has_locrests,
        seed, label,
        energy, energies, lablen
      );
    }
    if (iscore == 1) continue;
    if (iscore == 2) break;
    print_struc_(
     seed,
     label, 
     energy,
     energies,
     nlig,
     ens,
     phi,
     ssi,
     rot,
     xa,
     ya,
     za,
     locrests,     
     morph,
     nhm,
     dlig,
     has_locrests,     
     lablen
    );
  }    
  extern char *shmlinks[100];
  extern int shmlinkcount;
  for (n = 0; n < shmlinkcount; n++) shm_unlink(shmlinks[n]);
  return 0;
}
