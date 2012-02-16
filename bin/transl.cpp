#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: $path/trans <PDB file> <x> <y> <z>\n");
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
  if (argc < 5) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  if (argc > 5) {
    fprintf(stderr, "Too many arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }  
  Coor *x; int natoms; 
  char **pdbstrings;  bool *pdblayout; int linecounter;
  FILE *fil = fopen(argv[1], "r");
  double dx = atof(argv[2]);
  double dy = atof(argv[3]);
  double dz = atof(argv[4]);  
  read_pdb2(fil,x,pdbstrings,pdblayout,natoms,linecounter);
  
  for (int n = 0; n < natoms; n++) {
    x[n][0]+=dx;x[n][1]+=dy;x[n][2]+=dz;
  }
  
  write_pdb2(stdout,x,pdbstrings,pdblayout,natoms,linecounter);
}
