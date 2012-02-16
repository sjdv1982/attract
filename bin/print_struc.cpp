#include "max.h"
#include <cstdio>
#include <cstring>

int nstruc = 0;

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
 const double *morph,
 const int *nhm,
 const modes2 &dlig, 
 int len_label
) 
{
  nstruc += 1;
  printf("#%d\n", nstruc);
  printf("### SEED %d\n", seed);
  if (len_label > 1) {
    char label[1000];
    memcpy(label, label0, len_label);
    label[len_label] = 0;
    printf("%s",label);
  }
  if (energies != NULL) {
    printf("## Energy: %.3f\n", energy);
    printf("##   %.3f %.3f %.3f %.3f %.3f %.3f\n", 
     energies[0],energies[1],energies[2],
     energies[3],energies[4],energies[5]);
  }
  for (int i = 0; i < nlig; i++) {
    printf("   ");   
    if (morph[i] >= 0) {
      printf("%.4f ", morph[i]);
    }    
    else if (ens[i] > 0) {
      printf("%d ", ens[i]);
    }
    printf("%.8f %.8f %.8f %.4f %.4f %.4f", 
      phi[i], ssi[i], rot[i], 
      xa[i],  ya[i], za[i]      
    );
    for (int ii = 0; ii < nhm[i]; ii++) {
      printf(" %.4f", dlig[i][ii]);
    }
    printf("\n"); 
  }
  //fflush(stdout);
}

extern "C" void print_struc2_(
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
 const double *morph,
 const int *nhm,
 const modes2 &dlig, 
 const int &len_label
) {

 print_struc_(
  seed,
  label0, 
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
  morph,
  nhm,
  dlig, 
  len_label
 );
}
