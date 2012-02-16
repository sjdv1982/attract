#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "max.h"
#include <cmath>

void usage() {
  fprintf(stderr, "usage: $path/deflate <PDB file> <step size> <ensemble size>\n");
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

Coor *xx[10000];
static Coor grad[TOTMAXATOM];

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

double rmsd_center(Coor *x, int natoms, Coor &center) {
  double sd = 0;
  for (int n = 0; n < natoms; n++) {
    Coor &cx = x[n];
    double dx = cx[0]-center[0];
    double dy = cx[1]-center[1];
    double dz = cx[2]-center[2];
    double dsq = dx*dx+dy*dy+dz*dz;
    sd += dsq;
  }
  double rmsd = sqrt(sd/natoms);
  return rmsd;
}

void deflate(Coor *x, Coor *x2, int natoms, double step) {
  memset(grad, 0, natoms*sizeof(Coor));
  for (int n = 0; n < natoms; n++) {
    Coor &cx1 = x[n];
    Coor &g1 = grad[n];
    for (int nn = n+1; nn < natoms; nn++) {
      Coor &cx2 = x[nn];    
      Coor &g2 = grad[nn];
      double dx = cx1[0] - cx2[0];
      double dy = cx1[1] - cx2[1];
      double dz = cx1[2] - cx2[2];
      double ddsq = 1/(dx*dx+dy*dy+dz*dz);
      if (ddsq > 2) ddsq = 2;
      dx *= ddsq; dy *= ddsq; dz *= ddsq;
      g1[0] -= dx; g2[0] += dx;
      g1[1] -= dy; g2[1] += dy;
      g1[2] -= dz; g2[2] += dz;
    }
  }
  double totgrad = 0;
  for (int n = 0; n < natoms; n++) {
    Coor &g = grad[n];
    totgrad += sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);    
  }
  totgrad /= natoms;
  double gradfac = step/totgrad;
  for (int n = 0; n < natoms; n++) {
    Coor &cx = x[n];
    Coor &cx2 = x2[n];
    Coor &g = grad[n];
    cx2[0] = cx[0] + g[0]*gradfac;
    cx2[1] = cx[1] + g[1]*gradfac;
    cx2[2] = cx[2] + g[2]*gradfac;
  }
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Too few arguments\n"); usage();
  }
  if (argc > 4) {
    fprintf(stderr, "Too many arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }  
  
  double step = atof(argv[2]);
  int ens = atoi(argv[3]);
  
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
    
  Coor *curr = x;
  int count = 0;  
  while (1) {
    double rmsd = rmsd_center(curr, natoms, center);
    if (rmsd < 2*step) break;
    fprintf(stderr, "RMSD: %.3f\n", rmsd);
    xx[count] = new Coor[natoms];
    deflate(curr, xx[count], natoms, step);
    curr = xx[count];
    count += 1;
  }
  printf("MODEL %d\n", count+1);
  write_pdb2(stdout,x,pdbstrings,pdblayout,natoms,linecounter);
  printf("ENDMDL\n");

  for (int n = 0; n < ens-1; n++) {  
    int pos = int(double(n+1)/ens*count+0.5-1);
    Coor *curr = xx[pos];
    printf("MODEL %d\n", n+1);
    write_pdb2(stdout,curr,pdbstrings,pdblayout,natoms,linecounter);
    printf("ENDMDL\n");
  }
  
}
