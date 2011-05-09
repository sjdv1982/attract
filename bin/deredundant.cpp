//Converts DOF so that receptor rotations and translations are zero

//usage: ./fix_receptor structures.dat <number of ligands>


#include "max.h"
#include <cmath>
const double pi = 4.0f * atan(1.0f);

const double radgyr = 30.0;
const double lim = 0.05;

extern "C" void read_dofs_(const char *f_, dof &phi, dof &ssi, dof &rot, dof &xa, dof &ya, dof &za, modes &dlig, const int &nlig, int *nhm, int (&seed)[MAXSTRUC], char *label[MAXSTRUC], int &nstruc, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);

/* DOFs */
static double phi[MAXSTRUC][MAXLIG];
static double ssi[MAXSTRUC][MAXLIG];
static double rot[MAXSTRUC][MAXLIG];
static double xa[MAXSTRUC][MAXLIG];
static double ya[MAXSTRUC][MAXLIG];
static double za[MAXSTRUC][MAXLIG];
static double dlig[MAXSTRUC][MAXLIG][MAXMODE];
static double rotmat[MAXSTRUC][MAXLIG][9];

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/deredundant structures.dat <number of ligands>\n");
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
  int n, i;
  if (argc != 3) {
    fprintf(stderr, "Wrong number of arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }
  int nlig = atoi(argv[2]);
  int nstruc;

  //read DOFs and set pivots
  double pivot[3][MAXLIG];
  memset(pivot,0,sizeof(pivot));
  int nhm[MAXLIG];
  for (int n = 0; n < MAXLIG; n++) nhm[n] = MAXMODE;
  int seed[MAXSTRUC];
  char *label[MAXSTRUC];
  int auto_pivot, centered_receptor, centered_ligands;    
  read_dofs_(argv[1], phi, ssi, rot, xa, ya, za, dlig, nlig, nhm, seed, label, nstruc, pivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  
  if (centered_receptor != centered_ligands) { 
    fprintf(stderr, "Receptor and ligands must be both centered or both uncentered\n");
    exit(1);
  }  
  if (!(centered_ligands) && auto_pivot) {
    fprintf(stderr, "With uncentered ligands, pivots must be supplied\n");
    exit(1);
  }

  //determine number of modes
  for (i = 0; i < nlig; i++) {
    for (int m = 0; m < MAXMODE; m++) {
      bool zeromode = 1;
      for (n = 0; n < nstruc; n++) {
        if (fabs(dlig[n][i][m]) > 0.001) {
	  zeromode = 0;
	  break;
	}
      }
      if (zeromode) {
        nhm[i] = m;
	break;
      }
    }
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
  int nr[MAXSTRUC];
  int nonredundant = 0;
  for (n = 0; n < nstruc; n++) {    
    for (i = 0; i < nlig; i++) {
      euler2rotmat_(phi[n][i],ssi[n][i],rot[n][i],rotmat[n][i]);
    }
  }
  for (n = 0; n < nstruc; n++) {    

    bool unique = 1;
    for (int nn = 0; nn < nonredundant; nn++) {  
      double rmsd = 0;
      bool different = 0;
      for (i = 1; i < nlig; i++) { //assume fixed receptor...
        int n2 = nr[nn];
        for (int j=0;j<3;j++) {
          double *ra1 = &rotmat[n][i][3*j];
	  double *ra2 = &rotmat[n2][i][3*j];
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
        double dt[3] = {xa[n][i]-xa[n2][i],ya[n][i]-ya[n2][i],za[n][i]-za[n2][i]};
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
  
    nr[nonredundant] = n;
    nonredundant++;
    printf("#%d\n",nonredundant);
    printf("##%d => deredundant\n",n+1);
    if (label[n] != NULL) printf("%s",label[n]);
    
    
    for (i = 0; i < nlig; i++) {
      printf("   %.8f %.8f %.8f %.4f %.4f %.4f", 
        phi[n][i], ssi[n][i], rot[n][i], 
        xa[n][i],  ya[n][i], za[n][i]      
      );    
      for (int ii = 0; ii < nhm[i]; ii++) {
	printf(" %.4f", dlig[n][i][ii]);
      }
      printf("\n");     
    }
  }
}
