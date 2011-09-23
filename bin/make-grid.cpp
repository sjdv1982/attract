#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "grid.h"

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

const int gridextension = 32;

extern int cartstate_new(int argc, char *argv[],bool single=0);
extern "C" void cartstate_pivot_auto_(const int &handle);

extern "C" void cartstate_translate_atomtypes_(const int &handle);
extern "C" void cartstate_get_haspar_(const int &handle,iParameters *&haspar);

extern void get_shm_name(int shm_id, char *shm_name);
extern bool exists(const char *f);

int new_shm_id() {
  int ret = 1;
  char shm_name[100];
  while (1) {
    get_shm_name(ret, shm_name);
    int result = shm_open(shm_name, (O_CREAT|O_EXCL), (S_IREAD | S_IWRITE));
    if (result != -1) {
      return ret;
    }
    ret++;
  }
}

int main(int argc, char*argv[]) {
  if (argc < 7) {
    fprintf(stderr,"Too few arguments\n");
    fprintf(stderr, "Please provide PDB file, SITUS map of the interior, ATTRACT parameter file, plateau distance, neighbour distance, output file name  [--alphabet <alphabet file>] [--shm] [--cdie] [--epsilon=15] [calc-potentials=1]\n");
    return -1;
  }
  bool use_shm = 0;
  bool cdie = 0;
  float epsilon = 15;
  
  char *alphabetfile = NULL;
  for (int n = 1; n < argc; n++) {
    if (n < argc -1) {
      if (!strcmp(argv[n],"--epsilon")) {
        epsilon = atof(argv[n+1]);
        char **argv2 = new char *[argc-1];
        if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
        if (n+2 < argc) memcpy(argv2+n,argv+n+2,(argc-n-2)*sizeof(char*));
        argv = argv2;
        argc -= 2;
        n -= 1;
        continue;
      }
      if (!strcmp(argv[n],"--alphabet")) {
        alphabetfile = argv[n+1];
        if (!exists(alphabetfile)) {
          fprintf(stderr, "Alphabet file %s does not exist", alphabetfile); 
          exit(1);
        }
        char **argv2 = new char *[argc-1];
        if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
        if (n+2 < argc) memcpy(argv2+n,argv+n+2,(argc-n-2)*sizeof(char*));
        argv = argv2;
        argc -= 2;
        n -= 1;
        continue;
      }
    }
    if (!strcmp(argv[n],"--shm")) {
      use_shm = 1;
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+1 < argc) memcpy(argv2+n,argv+n+1,(argc-n-1)*sizeof(char*));
      argv = argv2;
      argc -= 1;
      n -= 1;
      continue;
    }
    if (!strcmp(argv[n],"--cdie")) {
      cdie = 1;
      char **argv2 = new char *[argc-1];
      if (n > 0) memcpy(argv2, argv,n*sizeof(char*));
      if (n+1 < argc) memcpy(argv2+n,argv+n+1,(argc-n-1)*sizeof(char*));
      argv = argv2;
      argc -= 1;
      n -= 1;
      continue;
    }
    
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

  bool alphabet[MAXATOMTYPES];
  memset(alphabet, 0, MAXATOMTYPES*sizeof(bool));  
  if (alphabetfile != NULL) {
    //custom alphabet:
    FILE *f = fopen(alphabetfile, "r");
    char buf[2000];
    while (!feof(f)) {
      if (!fgets(buf, 2000, f)) break;
      int residue;
      sscanf(buf, "%d", &residue);
      if (residue <= 0 || residue > MAXATOMTYPES) {
        fprintf(stderr, "Reading error in alphabet file: %s\n", buf);
        exit(1);
      }
      alphabet[residue-1] = 1;
    }
  }
  else {
    //default alphabet: 
    for (int i = 0; i < MAXATOMTYPES; i++) { 
      for (int j = 0; j < MAXATOMTYPES; j++) {
        if (haspar[i][j]) {
          alphabet[i] = 1;
          break;
        }
      }
    }
    if (alphabet[MAXATOMTYPES-1]) {
      fprintf(stderr, "WARNING: your parameter file is full-size and you didn't supply an alphabet file\nGrid calculation will take very long!\n");
    }
  }
  
  bool calc_pot = 1;
  if (argc > 7) calc_pot = atoi(argv[7]);
  g.calculate(cartstatehandle, 0, argv[2], atof(argv[4]), atof(argv[5]), 
   gridextension, 0, alphabet, cdie, epsilon, calc_pot);
  if (use_shm) {
    g.shm_energrads = new_shm_id();
    g.shm_neighbours = new_shm_id();  
  }
  g.write(argv[6]);
}  



