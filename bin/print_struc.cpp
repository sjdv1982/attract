#include "max.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;


int nstruc = 0;
bool has_flushed = 0;
time_t last_flush;

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
 const int *nihm,
 const modes2 &dlig, 
 const int *has_locrests,
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
    printf("%.6f %.6f %.6f %.4f %.4f %.4f", 
      phi[i], ssi[i], rot[i], 
      xa[i],  ya[i], za[i]      
    );
    if (has_locrests[i]) {
      for (int ii = 0; ii < 3; ii++) {
        printf(" %.4f", locrests[i][ii]);
      }
    }
    for (int ii = 0; ii < nhm[i]; ii++) {
      printf(" %.4f", dlig[i][ii]);
    }
    for (int ii = nhm[i]; ii < nhm[i]+nihm[i]; ii++) {
          printf(" %.4f", dlig[i][ii]);
        }
    printf("\n"); 
  }
  time_t curr_time = time(NULL);
  if (!has_flushed || difftime(curr_time, last_flush) > 0.5) {    
    fflush(stdout);
    has_flushed = 1;
    last_flush = curr_time;
  }  
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
 const coors2 &locrests,
 const double *morph,
 const int *nhm,
 const int *nihm,
 const modes2 &dlig, 
 const int *has_locrests,
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
  locrests,
  morph,
  nhm,
  nihm,
  dlig, 
  has_locrests,
  len_label
 );
}
