#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/center <PDB file>\n");
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

typedef double Coor[3];

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

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  if (argc > 2) {
    fprintf(stderr, "Too many arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }  
  Coor *x; int natoms; 
  char **pdbstrings;  bool *pdblayout; int linecounter;
  FILE *fil = fopen(argv[1], "r");
  read_pdb2(fil,x,pdbstrings,pdblayout,natoms,linecounter);
  
  Coor center = {0,0,0};
  for (int n = 0; n < natoms; n++) {
    center[0] += x[n][0]; center[1] += x[n][1]; center[2] += x[n][2];
  }
  if (natoms) {
    center[0]/=natoms;center[1]/=natoms;center[2]/=natoms;
  }
  fprintf(stderr, "Center: %.3f %.3f %.3f\n", center[0],center[1],center[2]);
  for (int n = 0; n < natoms; n++) {
    x[n][0]-=center[0];x[n][1]-=center[1];x[n][2]-=center[2];
  }
  
  write_pdb2(stdout,x,pdbstrings,pdblayout,natoms,linecounter);
}
