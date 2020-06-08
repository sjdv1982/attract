#include "state.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
using namespace std;

bool is_empty_(char *buf) {
  for (int n = 0; n < strlen(buf)-1; n++) {
    if (buf[n] != 32) return 0;
  }
  return 1;
}

void parse_restraintfile(MiniState &ms, const char *restfile) {  
  
  map<string,int *> selections;  
  map<string,int> sel_natoms;
  static Restraint restraints[MAXRESTRAINTS];
  int nr_restraints = 0;
  if (ms.nr_restraints > 0) {
    memcpy(restraints, ms.restraints, ms.nr_restraints * sizeof(Restraint));
    nr_restraints = ms.nr_restraints;
    delete[] ms.restraints;
  }
  
  char chars[] = {32,10,13,0};    
  
  FILE *f = fopen(restfile, "r");  
  int mode = 0; //0 = selection section, 1 = restraints section
  int line = 0;
  while (!feof(f)) {
    line++;
    char buf[100000];
    char l[100000];
    char *fields[10000];
    if (!fgets(buf,100000,f) || (!strlen(buf)) || is_empty_(buf)) {
      if (mode == 1) break;
      mode = 1; continue;
    }
    if (buf[0] == '#') continue;
    strcpy(l,buf);
    char *field = strtok(buf,chars);
    int nf = 0;
    while (field != NULL) {
      if (nf == 10000) {
        fprintf(stderr, "Reading error in %s, line %d: too many values on a line\n", restfile, line);
        exit(1);
      }
      fields[nf] = field;
      field = strtok(NULL, chars);
      nf++;
    }
    
    if (mode == 0) { //new selection
      if (nf < 2) {
        fprintf(stderr, "Reading error in %s, line %d: too few fields\n", restfile, line);
        exit(1);      
      }
      string id(fields[0]);
      int natoms = atoi(fields[1]);
      if (natoms > MAXSELECTION) {
        fprintf(stderr, "Reading error in %s, line %d: selection is too large\n", restfile, line);      
	exit(1);
      }
      if (nf != natoms + 2) {
        fprintf(stderr, "Reading error in %s, line %d %d: wrong number of atoms %d in selection\n", restfile, line, nf, natoms);      
	exit(1);	
      }
      int *atoms = new int[natoms+1];
      atoms[natoms] = 0;
      for (int n = 0; n < natoms; n++) atoms[n] = atoi(fields[n+2]);
      selections[id] = atoms;
      sel_natoms[id] = natoms;
    }
    else  { //new restraint
      if (nr_restraints == MAXRESTRAINTS) {
        fprintf(stderr, "Reading error in %s, line %d: maximum number of restraints exceeded\n", restfile, line);
	exit(1);      
      }
      if (nf < 3) {
        fprintf(stderr, "Reading error in %s, line %d: too few fields\n", restfile, line);
	exit(1);      
      }
      if (nf > 9) {
        fprintf(stderr, "Reading error in %s, line %d: too many fields\n", restfile, line);
        exit(1);      
      }
      string id1(fields[0]);
      if (selections.find(id1) == selections.end()) {
        fprintf(stderr, "Reading error in %s, line %d: Unknown selection %s\n", restfile, line, id1.c_str());
	exit(1);      
      }            
      Restraint r;
      memset(&r, 0, sizeof(r));
      r.selection1 = selections[id1];
      r.s1 = sel_natoms[id1];
      int maxindex = 0;
      for (int n = 0; n < r.s1; n++) {
        if (r.selection1[n] > maxindex) maxindex = r.selection1[n];
      }

      if (strlen(fields[1]) == 1 && fields[1][0] == '7') {
        r.type = 7;
      }
      else {
        r.type = atoi(fields[2]);
        if (r.type > 10 || r.type == 7) {
          fprintf(stderr, "Reading error in %s, line %d: Restraint type %s not supported\n", restfile, line, fields[2]);
          exit(1);      
        }
        string id2(fields[1]);
        if (selections.find(id2) == selections.end()) {
          fprintf(stderr, "Reading error in %s, line %d: Unknown selection %s\n", restfile, line, id2.c_str());
          exit(1);      
        }
        r.selection2 = selections[id2];
        r.s2 = sel_natoms[id2];
        for (int n = 0; n < r.s2; n++) {
          if (r.selection2[n] > maxindex) maxindex = r.selection2[n];
        }
      }  
      r.maxindex = maxindex;

      if (r.type > 10) {
        fprintf(stderr, "Reading error in %s, line %d: Restraint type %s not supported\n", restfile, line, fields[2]);
	exit(1);      
      }

      bool wrong = 0;
      if (r.type == 1 && nf != 5) wrong = 1;
      if (r.type == 2 && nf != 7) wrong = 1;
      if (r.type == 3 && nf != 5) wrong = 1;
      if (r.type == 4 && nf != 5) wrong = 1;
      if (r.type == 5 && nf != 5) wrong = 1;
      if (r.type == 6 && nf != 6) wrong = 1;
      if (r.type == 7 && nf != 9) wrong = 1;
      if (r.type == 8 && nf != 7) wrong = 1;
      if (r.type == 9 && nf != 5) wrong = 1;
      if (r.type == 10 && nf != 5) wrong = 1;
      if (wrong) {
        fprintf(stderr, "Reading error in %s, line %d: Wrong number of parameters for restraint type %d\n", restfile, line, r.type);
	exit(1);      
      }
      if (r.type == 3) {
        if (r.s1 != 1 || r.s2 != 1) {
          fprintf(stderr, "Reading error in %s, line %d: Minimum-distance restraints require single-atom selections\n", restfile, line);
	  exit(1);      
        }
      }
      if (r.type == 4) {
        if (r.s1 != 1 || r.s2 != 1) {
          fprintf(stderr, "Reading error in %s, line %d: harmonic bond restraints require single-atom selections\n", restfile, line);
      	  exit(1);
        }
      }
      if (r.type == 5) {
        if (r.s1 != 1 || r.s2 != 1) {
          fprintf(stderr, "Reading error in %s, line %d: LJ restraints require single-atom selections\n", restfile, line);
          exit(1);
        }
      }
      if (r.type == 6) {
        if (r.s1 != 1 || r.s2 != 1) {
          fprintf(stderr, "Reading error in %s, line %d: step potential requires single-atom selections\n", restfile, line);
          exit(1);
        }
      }   
      if (r.type == 8) {
        if (r.s1 != 1 || r.s2 != 1) {
          fprintf(stderr, "Reading error in %s, line %d: bump potential requires single-atom selections\n", restfile, line);
          exit(1);
        }
      } 
      if (r.type == 9) {
        if (r.s1 != 3 ) {
          fprintf(stderr, "Reading error in %s, line %d: angle potential requires 3-atom selections\n", restfile, line);
          exit(1);
        }
      } 
      if (r.type == 10) {
        if (r.s1 != 4 ) {
          fprintf(stderr, "Reading error in %s, line %d: torsion potential requires 4-atom selections\n", restfile, line);
          exit(1);
        }
      } 
      if (r.type == 7) {
        r.par1 = atof(fields[2]); //mindist
        r.par2 = atof(fields[3]); //maxdist
        r.par3 = atof(fields[4]); //force constant        
        string position_type = fields[5];
        for (int n = 0; n < position_type.length(); n++) {
          char c = position_type[n];          
          switch (c) {
            case 'x': 
              r.position_type |= 1;
              break;
            case 'y': 
              r.position_type |= 2;
              break;
            case 'z':
              r.position_type |= 4;
              break;              
            default:
              fprintf(stderr, "Reading error in %s, line %d: Position restraint type characters must be [xyz]\n", restfile, line);
              exit(1);            
          }
        }
        r.par4 = atof(fields[6]);  //refe-x
        r.par5 = atof(fields[7]);  //refe-y
        r.par6 = atof(fields[8]);  //refe-z
      }
      else {
        if (nf > 3) r.par1 = atof(fields[3]);
        if (nf > 4) r.par2 = atof(fields[4]);
        if (nf > 5) r.par3 = atof(fields[5]);
        if (nf > 6) r.par4 = atof(fields[6]);
        if (nf > 7) r.par5 = atof(fields[7]);
      }
      restraints[nr_restraints] = r;
      nr_restraints++;
    }
  }  
  ms.restraints = new Restraint[nr_restraints+1];
  ms.nr_restraints = nr_restraints;
  memcpy(ms.restraints, restraints, nr_restraints*sizeof(Restraint));

}
