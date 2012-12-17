//Converts DOF so that receptor rotations and translations are zero

//usage: ./fix_receptor structures.dat <number of ligands> [--ens <ensemble size for each ligand>] [--ignorens]\n");

//if ignorens == False (default), then different ensemble structures are never redundant
//if ignorens == True, then the ensemble index is simply ignored


#include "max.h"
#include <cmath>
#include <cstdio>

const double radgyr = 30.0;
const double lim = 0.05;

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

extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, coors2 &locrests, dof2 &morph, modes2 &dlig, const int &nlig, const int *nhm, const int *nrens0, const int *morphing, const int *has_locrests, int &seed, char *&label, int f_len);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

/* DOFs */
static int nrens[MAXLIG];
static int nhm[MAXLIG];

static int cens[MAXLIG];
static double cphi[MAXLIG];
static double cssi[MAXLIG];
static double crot[MAXLIG];
static double cxa[MAXLIG];
static double cya[MAXLIG];
static double cza[MAXLIG];
static double cmorph[MAXLIG];
static double dlig[MAXLIG][MAXMODE];

static double crotmat[MAXLIG][9];

static int ens[MAXSTRUC][MAXLIG];
static double phi[MAXSTRUC][MAXLIG];
static double ssi[MAXSTRUC][MAXLIG];
static double rot[MAXSTRUC][MAXLIG];
static double xa[MAXSTRUC][MAXLIG];
static double ya[MAXSTRUC][MAXLIG];
static double za[MAXSTRUC][MAXLIG];
static double rotmat[MAXSTRUC][MAXLIG][9];

static int seed;
static char *label;

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/deredundant structures.dat <number of ligands> [--modes <mode file>] [--ens <ensemble size for each ligand>] [--ignorens]\n");
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
  
  bool ignorens = 0;
  
  for (int n = 0; n < MAXLIG; n++) {
    nhm[n] = 0;
    nrens[n] = 0;
  }

  coors2 clocrests;
  int has_locrests[MAXLIG];
  memset(has_locrests, 0, MAXLIG*sizeof(int));

  while (argc > 3) {
    if (!strcmp(argv[3],"--ignorens")) {
      ignorens = 1;
      memmove(argv+3, argv+4, sizeof(char*) * (argc-3));
      argc--;  
      continue;    
    }
    if (!strcmp(argv[3],"--modes")) {
      int count = 0;
      while (argc > 4) {
        memmove(argv+3, argv+4, sizeof(char*) * (argc-3));
        argc--;      
        if (!strncmp(argv[3],"--",2)) break;      
        nhm[count] = atoi(argv[3]);
        count++;
      }
      argc--;
      continue;          
    }
    if (!strcmp(argv[3],"--ens")) {
      int count = 0;
      while (argc > 4) {
        memmove(argv+3, argv+4, sizeof(char*) * (argc-3));
        argc--;      
        if (!strncmp(argv[3],"--",2)) break;      
        nrens[count] = atoi(argv[3]);
        count++;
      }
      argc--;
      continue;
    }
    if (argc > 4 && (!strcmp(argv[3],"--locrest"))) {
      int lig = atoi(argv[4]);
      if (lig <= 0 || lig > MAXLIG) {
        fprintf(stderr,"Ligand %d must be larger than 0\n", lig);
        usage();
      }
      has_locrests[lig-1] = 1;
      memmove(argv+3, argv+5, sizeof(char*) * (argc-4));
      argc -= 2;      
      continue;
    }
    fprintf(stderr, "Wrong number of arguments\n"); usage();
  }  

  if (argc != 3) {
    fprintf(stderr, "Wrong number of arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }
  
  int nlig = atoi(argv[2]);

  //read DOFs and set pivots
  double pivot[3][MAXLIG];
  memset(pivot,0,sizeof(pivot));
  
  int auto_pivot, centered_receptor, centered_ligands;    
  int line;
  FILE *fil = read_dof_init_(argv[1], nlig, line, pivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1])); 
  if (centered_receptor != centered_ligands) { 
    fprintf(stderr, "Receptor and ligands must be both centered or both uncentered\n");
    exit(1);
  }  
  if (!(centered_ligands) && auto_pivot) {
    fprintf(stderr, "With uncentered ligands, pivots must be supplied\n");
    exit(1);
  }

  if (auto_pivot) printf("#pivot auto\n");
  else {
    for (i = 0; i < nlig; i++) {
      printf("#pivot %d %.3f %.3f %.3f\n", 
	i+1, pivot[0][i], pivot[1][i], pivot[2][i]);
    }   	
  }   
  if (centered_receptor) printf("#centered receptor: true\n");
  else printf("#centered receptor: false\n");
  if (centered_ligands) printf("#centered ligands: true\n");
  else printf("#centered ligands: false\n");
  //main loop  
  int nonredundant = 0;
  int nstruc = 0;

  int morphing[MAXLIG];
  memset(morphing, 0, MAXLIG*sizeof(int));

  while (1) {

    int result = read_dof_(fil, line, nstruc, argv[1], cens, cphi, cssi, crot, 
     cxa, cya, cza, clocrests, 
     cmorph, dlig, nlig, nhm, nrens, morphing, has_locrests, 
     seed, label, strlen(argv[1])
    );
    if (result != 0) break;

    if ((fabs(cphi[0])>0.001)|| (fabs(cssi[0])>0.001) ||(fabs(crot[0])>0.001)||
        (fabs(cxa[0])>0.001) || (fabs(cya[0])>0.001) ||(fabs(cza[0])>0.001))
    {
      fprintf(stderr, "Structure %d: receptor not fixed\n", nstruc);
      exit(1);  
    }
    for (i = 0; i < nlig; i++) {
      euler2rotmat_(cphi[i],cssi[i],crot[i],crotmat[i]);
    }

    bool unique = 1;
    for (int nn = 0; nn < nonredundant; nn++) {  
      double rmsd = 0;
      bool different = 0;
      for (i = 0; i < nlig; i++) { 
        if (nrens[i] && (cens[i] != ens[nn][i]) && ignorens == 0) {
          different = 1;
          break;
        }
        if (i == 0) continue; //assume fixed receptor...
        
        for (int j=0;j<3;j++) {
          double *ra1 = &crotmat[i][3*j];
	  double *ra2 = &rotmat[nn][i][3*j];
	  double dr[3] = {ra1[0]-ra2[0],ra1[1]-ra2[1],ra1[2]-ra2[2]};
	  double dissq = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
	  rmsd += radgyr/(nlig-1) * dissq/3;
          if (rmsd > lim) {
  	    //printf("DIF1 %d %d %.3f\n", n,n2,rmsd);
            different = 1;
	    break;
          }	
	}
        if (different) break;	  
        double dt[3] = {
         cxa[i]-xa[nn][i],
         cya[i]-ya[nn][i],
         cza[i]-za[nn][i]
        };
	//printf("%.3f %.3f %.3f\n", dt[0],dt[1],dt[2]);
	double dissq = dt[0]*dt[0]+dt[1]*dt[1]+dt[2]*dt[2];
	rmsd += dissq/(nlig-1);
        if (rmsd > lim) {
	  //printf("DIF2 %d %d %.3f\n", n,n2,rmsd);
          different = 1;
	  break;
        }	
      }
      if (!different) {
        unique = 0;
	break;
      }
    }
    if (!unique) continue;
      
    if (nonredundant == MAXSTRUC) {
      fprintf(stderr, "Too many non-redundant structures. Increase MAXSTRUC and recompile\n");
      exit(1);
    }  
    memcpy(rotmat[nonredundant], crotmat, sizeof(crotmat));
    ens[nonredundant][0] = cens[0];
    for (i = 1; i < nlig; i++) {    
      ens[nonredundant][i] = cens[i];
      phi[nonredundant][i] = cphi[i];
      ssi[nonredundant][i] = cssi[i];
      rot[nonredundant][i] = crot[i];
      xa[nonredundant][i] = cxa[i];
      ya[nonredundant][i] = cya[i];
      za[nonredundant][i] = cza[i];
    }
    nonredundant++;
    
    char extralabel[1000];
    sprintf(extralabel, "### %d => deredundant\n", nstruc);
    double dummy = 0;
    int lablen = 0;
    if (label) {
      char label2[1000];
      strcpy(label2, label);
      strcat(label2, extralabel);
      label = label2;
      lablen = strlen(label);
    }
    else {
      label = extralabel;
      lablen = strlen(extralabel);
    }

    print_struc_(
     seed,
     label, 
     dummy,
     NULL,
     nlig,
     cens,
     cphi,
     cssi,
     crot,
     cxa,
     cya,
     cza,
     clocrests,
     cmorph,
     nhm,
     dlig,
     has_locrests,
     lablen
    );
  }
}
