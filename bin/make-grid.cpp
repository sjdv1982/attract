#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "grid.h"

const int gridextension = 32;

extern int cartstate_new(int argc, char *argv[],bool single=0);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_translate_atomtypes_(const int &handle);
extern "C" void cartstate_get_haspar_(const int &handle,iParameters *&haspar);

int main(int argc, char*argv[]) {
  if (argc < 7) {
    fprintf(stderr,"Too few arguments\n");
    fprintf(stderr, "Please provide PDB file, SITUS map of the interior, ATTRACT parameter file, plateau distance, neighbour distance, output file name\n");
    return -1;
  }
   //load the Cartesian parameters and get a handle to it
  int cartstatehandle;
  char *argv0[] = {argv[3], argv[1]};
  cartstatehandle = cartstate_new(2, argv0, 1);
  cartstate_pivot_auto_(cartstatehandle);
  Grid g;
  cartstate_translate_atomtypes_(cartstatehandle);
  iParameters *haspar0;
  cartstate_get_haspar_(cartstatehandle, haspar0);
  iParameters &haspar = *haspar0;

  //default alphabet: 
  bool alphabet[MAXATOMTYPES];
  memset(alphabet, 0, MAXATOMTYPES*sizeof(bool));  
  for (int i = 0; i < MAXATOMTYPES; i++) { 
    for (int j = 0; j < MAXATOMTYPES; j++) {
      if (haspar[i][j]) {
        alphabet[i] = 1;
        break;
      }
    }
  }

  //TODO: allow non-default alphabet

  g.calculate(cartstatehandle, 0, argv[2], atof(argv[4]), atof(argv[5]), gridextension, 0, alphabet);
  g.write(argv[6]);
}  
