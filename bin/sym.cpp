#include "max.h"
#include "state.h"
#include <cstdio>
#include <cstdlib>

#define MAXSYMHANDLES 10

static int nhandles = 0;
static int handles[MAXSYMHANDLES];

static int amino[MAXSYMHANDLES][MAXLIG][MAXRES];
static int namino[MAXSYMHANDLES][MAXLIG];

extern CartState &cartstate_get(int handle);

const double symfactor = 0.005;

int find_amino(int *aminolist, const int *iaci, int ncoor) {
  int namino = 0;
  for (int n = 0; n < ncoor; n++) {
    if (iaci[n] == 30) {
      aminolist[namino] = n;
      namino++;
    }
  } 
  return namino;
}

void test(int cartstatehandle) {
  if (nhandles == MAXSYMHANDLES) {
    fprintf(stderr, "Too many cartstates with symmetries\n");
    exit(1);
  }
  CartState &cs = cartstate_get(cartstatehandle);
  for (int n = 0; n < cs.nsym;n++) {
    int symtype = cs.symtypes[n];
    int ncoor;
    int firstlig;
    for (int nn = 0; nn < symtype; nn++) {
      int lig = cs.sym[n][nn]-1;
      int cncoor;
      int *iaci;
      int camino[MAXRES];


      if (lig == 0) {
        cncoor = cs.ieins[0];
        iaci = cs.iaci;
      }
      else {
        cncoor = cs.ieins[lig]-cs.ieins[lig-1];      
        iaci = &cs.iaci[cs.ieins[lig-1]];      
      }
      
      if (nn == 0) {
        firstlig = lig;
        ncoor = cncoor;
        namino[nhandles][n] = find_amino(amino[nhandles][n], iaci, ncoor);      
        if (namino[nhandles][n] == 0) {
          fprintf(stderr, "Ligand %d is in a symmetry group, but doesn't contain any residues (0 amino atoms)\n", lig+1);
          exit(1);
        }
      }
      else {
        if (cncoor != ncoor) {
          fprintf(stderr, "Ligands %d and %d are in the same symmetry group, but they don't have the same number of atoms: %d vs %d\n", firstlig+1,lig+1, ncoor, cncoor);
          exit(1);
        }        
        int cnamino = find_amino(camino, iaci, cncoor);
        if (cnamino != namino[nhandles][n]) {
          fprintf(stderr, "Ligands %d and %d are in the same symmetry group, but they don't have the same number of residues: %d vs %d\n", firstlig+1,lig+1, namino[nhandles][n], cnamino);
          exit(1);
        }        
        for (int nnn = 0; nnn < namino[nhandles][n]; nnn++) {
          if (camino[nnn] != amino[nhandles][n][nnn]) {
            fprintf(stderr, "Ligands %d and %d are in the same symmetry group, but they don't have their atoms in the same order\n", firstlig+1,lig+1);
            exit(1);
          }
        }        
      }  //end if
    } //next ligand   
  } // next symmetry group
  handles[nhandles] = cartstatehandle;
  nhandles++;
}

double symrest(int cyclesize, Coor **cx1, Coor **cf1, Coor **cx2, Coor **cf2) {
  double energy = 0;
  double avg = 0;
  for (int n = 0; n < cyclesize; n++) {
    Coor &x1 = *(cx1[n]);
    Coor &x2 = *(cx2[n]);
    double avgx = x1[0]-x2[0]; 
    double avgy = x1[1]-x2[1]; 
    double avgz = x1[2]-x2[2];     
    double avgsq =  avgx*avgx+avgy*avgy+avgz*avgz;  
    double avg0 = sqrt(avgsq);
    avg += avg0;
  }
  avg /= cyclesize; 
  

  for (int n = 0; n < cyclesize; n++) {
    Coor &x1 = *(cx1[n]);
    Coor &x2 = *(cx2[n]);
    Coor &f1 = *(cf1[n]);
    Coor &f2 = *(cf2[n]);
    double dx = x1[0] - x2[0];
    double dy = x1[1] - x2[1];
    double dz = x1[2] - x2[2];
    double dsq = dx*dx+dy*dy+dz*dz; 
    double d = sqrt(dsq);
    //printf("AVG %.3f D %.3f\n", avg, d);
    double dd = d - avg;
    double ene = symfactor * dd * dd;
    double fac = 2 * symfactor * dd;
    energy += ene;
    double fx = fac * dx/d;
    double fy = fac * dy/d;
    double fz = fac * dz/d;
    f1[0] -= fx; f1[1] -= fy; f1[2] -= fz;
    f2[0] += fx; f2[1] += fy; f2[2] += fz;    
  }
  return energy;
}

extern "C" void sym_(const int &cartstatehandle, const int &iab, double &energy) {

  CartState &cs = cartstate_get(cartstatehandle);
  if (cs.nsym == 0) return;
  
  int handleindex = -1;
  for (int n = 0; n < nhandles; n++) {
    if (handles[n] == cartstatehandle) {
      handleindex = n;
      break;
    }
  }
  if (handleindex == -1) {
    handleindex = nhandles;
    test(cartstatehandle);
  }
  //printf("DONE TESTING\n");
  
  Coor *x = (Coor *) &cs.x[0];
  Coor *f = (Coor *) &cs.f[0];

  for (int n = 0; n < cs.nsym; n++) {
    int symtype = cs.symtypes[n];  
    //printf("SYMTYPE %d\n", symtype);
    int na = namino[handleindex][n];
    for (int i = 0; i < na; i++) {
      int i2 = na-i-1;
      int aminoindex = amino[handleindex][n][i];
      int aminoindex2 = amino[handleindex][n][i2];
      Coor *cx1[MAXLIG];
      Coor *cx2[MAXLIG];
      Coor *cf1[MAXLIG];
      Coor *cf2[MAXLIG];
      int maxstep = ceil(symtype/2);
      for (int ligstep = 1; ligstep <= maxstep; ligstep++) {
        for (int nn = 0; nn < symtype; nn++) {
          int lig = cs.sym[n][nn]-1;
          int nextnn = nn + ligstep;
          if (nextnn >= symtype) nextnn -= symtype;
          int nextlig = cs.sym[n][nextnn]-1;

          int start = 0;
          if (lig > 0) start = cs.ieins[lig-1];
          int start2 = 0;
          if (nextlig > 0) start2 = cs.ieins[nextlig-1];
          //printf("%d START %d AMINOINDEX %d START2 %d AMINOINDEX2 %d\n", nn, start, aminoindex, start2, aminoindex2);
          cx1[nn] = &x[start+aminoindex];
          cf1[nn] = &f[start+aminoindex];
          cx2[nn] = &x[start2+aminoindex2];
          cf2[nn] = &f[start2+aminoindex2];        
        }
        double newenergy = symrest(symtype, cx1, cf1, cx2, cf2);
        energy += newenergy;
        //if (n == 1) break;
      }
    }
  }
  
  return;
}
