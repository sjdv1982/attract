#include "max.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory>

typedef double Coor[3];

void read_pdb2(
  FILE *fil, Coor *&x, 
  char **&pdbstrings, bool *&pdblayout,
  int &coorcounter, int &linecounter
) 
{
  Coor *x0 = new Coor[TOTMAXATOM];
  char **pdbstrings0 = new char*[3*TOTMAXATOM];
  bool *pdblayout0 = new bool[2*TOTMAXATOM];
  char buf[2000];
  
  linecounter = 0;
  coorcounter = 0;
  int stringcounter = 0;
  
  while (!feof(fil)) {
    char code[10];
    if (!fgets(buf, 2000, fil)) break;
    sscanf(buf, "%10s ", code);
    
    if (!strncmp(code,"ATOM", 4) || !strncmp(code,"HETATM", 8)) {
      Coor &c = x0[coorcounter];
      c[0] = atof(buf+30);
      c[1] = atof(buf+38);
      c[2] = atof(buf+46);
      coorcounter++;
      
      char *&s1 = pdbstrings0[stringcounter];      
      s1 = new char[31];
      memcpy(s1, buf, 30);
      s1[30] = 0;
      
      char *&s2 = pdbstrings0[stringcounter+1];
      char *end = buf+54;
      s2 = new char[strlen(end)+1];
      strcpy(s2, end);

      stringcounter += 2;
      
      pdblayout0[linecounter] = 1;
    }
    else {
      char *&s1 = pdbstrings0[stringcounter];
      s1 = new char[strlen(buf)+1];
      strcpy(s1, buf);

      stringcounter += 1;
    
      pdblayout0[linecounter] = 0;
    }
    linecounter++;
  }
  if (linecounter) {
    pdblayout = new bool[linecounter];
    memcpy(pdblayout, pdblayout0, linecounter*sizeof(bool));
  }
  if (coorcounter) {
    x = new Coor[coorcounter];
    memcpy(x,x0,coorcounter*sizeof(Coor));
  }
  if (stringcounter) {
    pdbstrings = new char*[stringcounter];
    memcpy(pdbstrings, pdbstrings0,stringcounter*sizeof(char *));
  }
  fclose(fil);
  delete [] x0; delete [] pdbstrings0; delete [] pdblayout0;
}

void write_pdb2(
  FILE *fil, const Coor *x, 
  char **pdbstrings, const bool *pdblayout,
  int coorcounter, int linecounter
) 
{
  int coorcounter0 = 0;
  for (int n = 0; n < linecounter; n++) {
    coorcounter0 += int(pdblayout[n]);
  }
  if (coorcounter0 != coorcounter) {
    fprintf(stderr, "write_pdb2: Wrong number of coordinates: %d given, %d computed\n", coorcounter, coorcounter0);
  }
    
  int currcoor = 0, currstring = 0;
  for (int n = 0; n < linecounter; n++) {
    if (pdblayout[n]) {
      fprintf(fil, "%30.30s%8.3f%8.3f%8.3f%s",
        pdbstrings[currstring], 
	x[currcoor][0],x[currcoor][1],x[currcoor][2],
	pdbstrings[currstring+1]
      );
      currcoor += 1;
      currstring += 2;
    }
    else {
      fprintf(fil, "%s", pdbstrings[currstring]);
      currstring += 1;
    }    
  }
}
