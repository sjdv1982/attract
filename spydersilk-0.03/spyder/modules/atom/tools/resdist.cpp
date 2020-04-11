// Copyright 2007, Sjoerd de Vries
// This file is part of the Spyder module: "atom" 
// For licensing information, see LICENSE.txt 

//DIFFERS from resdist in WHISCY because it does not use a conversion table

#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

struct Coor3f {
  float x;
  float y;
  float z;
};

struct Residue {
  int nr;
  vector<Coor3f> coor;
};

vector<Residue> res;

int main(int argc, char *argv[]) {
  char buf[2000];

  if (argc < 2) {
    fprintf(stderr,"ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: resdist <pdb file>\n");
    return 1;
  }

  char *filename = argv[1];

  FILE *fil = fopen(filename, "r");
  if (fil == NULL) {
    fprintf(stderr, "ERROR: PDB file %s does not exist\n", filename);
    return 1;
  }
  int currnr = -1;
  while (!feof(fil)) {
    char code[10];
    char atom[5];
    if (!fgets(buf, 2000, fil)) break;
    sscanf(buf, "%s %*d %s", code, atom);
    if (!strncmp(code,"ATOM", 4) && atom[0] != 'H') {
      int nr = atoi(buf + 22);
      if (nr != currnr) {
        Residue r;
        r.nr = nr;
	res.push_back(r);
	currnr = r.nr;
      }
      Residue &rcurr = res[res.size() -1];
      Coor3f ccoor;
      ccoor.x = atof(buf+27);
      ccoor.y = atof(buf+38);
      ccoor.z = atof(buf+46);
      rcurr.coor.push_back(ccoor);
    }
  }

  if (!res.size()) {fprintf(stderr, "ERROR: PDB file %s contains no residues\n", filename); return 1;}

  double mindissq0 = 10e10;

  for (int n = 0; n < res.size(); n++) {
    vector<Coor3f> &c1 = res[n].coor;
    for (int nn = n + 1; nn < res.size(); nn++) {
      vector<Coor3f> &c2 = res[nn].coor;
      double mindissq = mindissq0;
      for (int i = 0; i < res[n].coor.size(); i++) {
        for (int ii = 0; ii < res[nn].coor.size(); ii++) {
	  double currdissq =
	    (c1[i].x - c2[ii].x) * (c1[i].x - c2[ii].x) +
	    (c1[i].y - c2[ii].y) * (c1[i].y - c2[ii].y) +
	    (c1[i].z - c2[ii].z) * (c1[i].z - c2[ii].z)
	   ;
	   if (currdissq < mindissq) mindissq = currdissq;
        }
      }
      printf ("%d %d %f\n", res[n].nr, res[nn].nr, sqrt(mindissq));
    }
  }
}
