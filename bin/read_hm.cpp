/*
  read_hmm.cpp
  Author: Sjoerd de Vries

  Parses harmonic mode file
    
  Meant to be called from Fortran by ATTRACT
  The harmonic modes are automatically normalized
    
  new variables:
    nhm: the number of harmonic modes for each ligand
  
  changed variables:
    eigl: if mode is multi, this must be a THREE-DIMENSIONAL array allocated in Fortran:
    Fortran dimensions: maxlig x maxmode x max3atom (3 times the max number of atoms)
    C dimensions: max3atom x  maxmode x maxlig
    if mode is not multi, it must be a maxmode x max3atom array 

  File format:
  <number of mode> <mode amplitude>
  <any number of lines with any number of values per line>
  <number of next mode> <mode amplitude>
  <any number of lines with any number of values per line>
  ...
  it is assumed that the mode numbers are 1,2,3,..., i, 1,2,3,..., j, 1,2,3, ...
    as soon it goes back to 1, this is interpreted as the start of a new ligand
  in order to specify a ligand as having no modes, add a line just containing 0 for that ligand

*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "max.h"

const int MAX3ATOM = 3 * MAXATOM;
typedef double (*eigptr_multi)[MAXLIG][MAXMODE][MAX3ATOM];
typedef double (&eigref_multi)[MAXLIG][MAXMODE][MAX3ATOM];
typedef double (*eigptr_mono)[MAXMODE][MAX3ATOM];
typedef double (&eigref_mono)[MAXMODE][MAX3ATOM];
typedef double (*index_eigptr_multi)[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE];
typedef double (&index_eigref_multi)[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE];
typedef double (*index_eigptr_mono)[MAXINDEXMODE][MAXLENINDEXMODE];
typedef double (&index_eigref_mono)[MAXINDEXMODE][MAXLENINDEXMODE];
typedef int (*index_eigptr_multi2)[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE];
typedef int (&index_eigref_multi2)[MAXLIG][MAXINDEXMODE][MAXLENINDEXMODE];
typedef int (*index_eigptr_mono2)[MAXINDEXMODE][MAXLENINDEXMODE];
typedef int (&index_eigref_mono2)[MAXINDEXMODE][MAXLENINDEXMODE];


inline void check_hm(const char *hmfile, const int *natom, int line, int currlig, int currmode, int currpos, double *eigl, int multi) {
  if (natom[currlig-1] != -1 && currpos != 3 * natom[currlig-1]) {
    fprintf(stderr, "Reading error in %s, line %d: Ligand %d, mode %d: read %d values, expected %d values\n", hmfile, line, currlig, currmode, currpos,  3 * natom[currlig-1]);
    exit(1);
  }
  //normalize
  eigref_multi eigl_multi = *((eigptr_multi)(eigl));
  eigref_mono eigl_mono = *((eigptr_mono)(eigl));  

  double eigval = 0;
  int n;
  for (n = 0; n < currpos; n++) {
    double v;
    if (multi) v = eigl_multi[currlig-1][currmode-1][n];
    else v = eigl_mono[currmode-1][n];
    eigval += v * v;
  }
  eigval = sqrt(eigval);
  for (n = 0; n < currpos; n++) {
    double *v;
    if (multi) v = &eigl_multi[currlig-1][currmode-1][n];
    else v = &eigl_mono[currmode-1][n];
    *v /= eigval;
  }
}

extern "C" void read_hm_(const char *hmfile_, const char *hmword_, const int &nlig, const int *natom, int *nhm, double (&vall)[MAXLIG][MAXMODE], double *eigl, const int &multi, int hmfile_len, int hmword_len) {
  char hmfile[1000];
  memcpy(hmfile, hmfile_, hmfile_len);
  hmfile[hmfile_len] = 0;
  char hmword[1000];
  memcpy(hmword, hmword_, hmword_len);
  hmword[hmword_len] = 0;
  
  eigref_multi eigl_multi = *((eigptr_multi)(eigl));
  eigref_mono eigl_mono = *((eigptr_mono)(eigl));
  //printf("%s %d %d %d %d %d %d\n", hmfile, nlig, *natom, *nhm, multi, hmfile_len, hmword_len); 
  if (nlig > MAXLIG) {
    fprintf(stderr, "Error in read_hm: number of %ss is larger than %d\n",hmword, MAXLIG);
    exit(1);
  }  
  if (multi == 0 && nlig != 1) {
    fprintf(stderr, "Error in read_hm: if multi is false, nlig must be 1\n");
    exit(1);
  }
  FILE *fil = fopen(hmfile, "r");
  if (fil == NULL) {
    fprintf(stderr, "Error in read_hm: cannot open %s\n" , hmfile);
    exit(1);
  }
  char buf[100000];
  float fields[10000];
  int currlig = 1;
  int currmode = 0;
  int currpos = 0;
  int line = 0;
  bool reading = 0;
  while (!feof(fil)) {
    line++;
    if(!fgets(buf,100000,fil)) continue;
    char chars[] = {32,10,13,0};    
    char *field = strtok(buf,chars);
    int nf = 0;
    while (field != NULL) {
      if (nf == 10000) {
        fprintf(stderr, "Reading error in %s, line %d: too many values on a line\n", hmfile, line);
        exit(1);
      }
      fields[nf] = atof(field);
      field = strtok(NULL, chars);
      nf++;
    }
    if (nf == 1 && fields[0] == 0) {
      //End of this mode
      if (currmode > 0) {
        check_hm(hmfile, natom,line,currlig,currmode, currpos, eigl, multi);
        fprintf(stderr,"%d modes read for %s %d\n", currmode, hmword, currlig); 
        nhm[currlig-1] = currmode;
      }
      //We have reached the next ligand
      //No modes for this ligand      
      currlig++;           
      currmode = 0;
      currpos = 0;
      reading = 0;
      continue;
    }
    if (nf == 2 && fields[0] == currmode+1) {
      //End of this mode
      if (currmode > 0) {
        check_hm(hmfile, natom,line,currlig,currmode, currpos, eigl, multi);
      }
      //We have reached the next mode
      currmode++;
      vall[currlig-1][currmode-1] = fields[1] * fields[1];  //square the amplitude
      if (currmode > MAXMODE) {
        fprintf(stderr, "Reading error in %s, line %d: Cannot read more than %d modes for %s %d\n", hmfile, line, MAXMODE, hmword, currlig);        
	exit(1);
      }
      currpos=0;
      reading = 1;
      continue;
    }
    if (nf == 2 && fields[0] == 1) {
      //End of this mode
      if (currmode > 0) {
        check_hm(hmfile, natom,line,currlig,currmode, currpos, eigl, multi);
        fprintf(stderr,"%d modes read for %s %d\n", currmode, hmword, currlig); 
	nhm[currlig-1] = currmode;
      }
      
      //We have reached the next ligand
      currlig++;
      if (currlig > nlig) {
        fprintf(stderr, "Reading error in %s, line %d: Cannot read more than %d %ss\n", hmfile, line, nlig,hmword);        
	exit(1);
      }

      currmode = 1;
      currpos = 0;
      reading = 1;
      continue;
    }  
    
    //reading data
    if (!reading) {
      fprintf(stderr, "Reading error in %s, line %d: Not expecting values here\n", hmfile, line);
    }
    if (currpos + nf >= MAX3ATOM) {
      fprintf(stderr, "Reading error in %s, line %d: %s %d, mode %d: More than %d values specified\n", hmfile, line, hmword, currlig, currmode, MAX3ATOM);
    }
    for (int n = 0; n < nf; n++) {
      if (multi) {
        eigl_multi[currlig-1][currmode-1][currpos+n] = fields[n]; 
      }
      else {
        eigl_mono[currmode-1][currpos+n] = fields[n]; 
      }
    }
    currpos += nf;    
  }
  
  if (currmode > 0) {
    check_hm(hmfile, natom,line,currlig,currmode, currpos, eigl, multi);
  }
  fprintf(stderr,"%d modes read for %s %d\n", currmode, hmword, currlig); 
  nhm[currlig-1] = currmode;
  if (currlig != nlig) {
    fprintf(stderr, "Reading error in %s, line %d: Read %d %ss, expected %d ligands\n", hmfile, line, currlig, hmword, nlig);        
    exit(1);
  }
  
  
}

extern "C" void read_indexmode_(const char *hmfile_, const char *hmword_, const int &nlig, int *nhm, int *eigl, double *val_eigl, const int &multi, int hmfile_len, int hmword_len) {
  char hmfile[1000];
  memcpy(hmfile, hmfile_, hmfile_len);
  hmfile[hmfile_len] = 0;
  char hmword[1000];
  memcpy(hmword, hmword_, hmword_len);
  hmword[hmword_len] = 0;
  int count = 0;

  index_eigref_multi2 eigl_multi = *((index_eigptr_multi2)(eigl));
  index_eigref_mono2 eigl_mono = *((index_eigptr_mono2)(eigl));
  index_eigref_multi val_eigl_multi = *((index_eigptr_multi)(val_eigl));
  index_eigref_mono val_eigl_mono = *((index_eigptr_mono)(val_eigl));
  //printf("%s %d %d %d %d %d %d\n", hmfile, nlig, *natom, *nhm, multi, hmfile_len, hmword_len);
  if (nlig > MAXLIG) {
    fprintf(stderr, "Error in read_indexmode: number of %ss is larger than %d\n",hmword, MAXLIG);
    exit(1);
  }
  if (multi == 0 && nlig != 1) {
    fprintf(stderr, "Error in read_indexmode: if multi is false, nlig must be 1\n");
    exit(1);
  }
  FILE *fil = fopen(hmfile, "r");
  if (fil == NULL) {
    fprintf(stderr, "Error in read_indexmode: cannot open %s\n" , hmfile);
    exit(1);
  }
  char buf[100000];
  int fields[10000];
  int currlig = 0;
  int currmode = 0;
  int currpos = 0;
  int line = 0;
  bool reading = 0;
  while (!feof(fil)) {
    line++;
    if(!fgets(buf,100000,fil)) continue;
    char chars[] = {32,10,13,0};
    char *field = strtok(buf,chars);
    int nf = 0;
    while (field != NULL) {
      if (nf == 10000) {
        fprintf(stderr, "Reading error in %s, line %d: too many values on a line\n", hmfile, line);
        exit(1);
      }
      fields[nf] = atoi(field);
      field = strtok(NULL, chars);
      nf++;
    }
    if (nf == 1 && fields[0] == 0) {
      //End of this mode
      if (currmode > 0) {
        fprintf(stderr,"%d index modes read for %s %d\n", currmode, hmword, currlig);
        nhm[currlig-1] = currmode;
      }
      //We have reached the next ligand
      //No modes for this ligand
      currlig++;
      currmode = 0;
      currpos = 0;
      reading = 0;
      continue;
    }
    if (nf == 1 && fields[0] == -1) {
      //End of this mode
      if (currmode > 0) {
        fprintf(stderr,"%d index modes read for %s %d\n", currmode, hmword, currlig);
	nhm[currlig-1] = currmode;
      }

      //We have reached the next ligand
      currlig++;
      if (currlig > nlig) {
        fprintf(stderr, "Reading error in %s, line %d: Cannot read more than %d %ss\n", hmfile, line, nlig,hmword);
	exit(1);
      }

      currmode = 0;
      currpos = 0;
      reading = 1;
      count = 0;
      // Initialize index modes
      for (int c=0; c< MAXLENINDEXMODE; c++){
	for (int c2=0; c2< MAXINDEXMODE; c2++){
	  if (multi){
		  eigl_multi[currlig-1][c2][c] = -1;
		  val_eigl_multi[currlig-1][c2][c] = 0;
		  }
	  else{
		  eigl_mono[c2][c] = -1;
		  val_eigl_mono[c2][c] = 0;
		  }
	}
      }
      
      continue;
    }

    //reading data
    count = 0;
    if (!reading) {
      fprintf(stderr, "Reading error in %s, line %d: Not expecting values here\n", hmfile, line);
    }
    if (currpos + nf >= MAX3ATOM) {
      fprintf(stderr, "Reading error in %s, line %d: %s %d, mode %d: More than %d values specified\n", hmfile, line, hmword, currlig, currmode, MAX3ATOM);
    }
    for (int n = 0; n < nf; n++) {
      if (multi) {
    	  if ( fields[n] != 0 ){
    		  if ( count < MAXLENINDEXMODE ){
    		  eigl_multi[currlig-1][currmode][count] = 3*(fields[n]-1);
    		  val_eigl_multi[currlig-1][currmode][count] = 1.0;
    		  eigl_multi[currlig-1][currmode+1][count] = 3*(fields[n]-1)+1;
    		  val_eigl_multi[currlig-1][currmode+1][count] = 1.0;
    		  eigl_multi[currlig-1][currmode+2][count] = 3*(fields[n]-1)+2;
    		  val_eigl_multi[currlig-1][currmode+2][count] = 1.0;
    		  //fprintf(stderr,"Nonzero %i %i %i\n", currlig-1, currmode, fields[n]);
    		  count ++;
    		  }
    		  else{
    			  fprintf(stderr, "Reading error in %s, line %d index mode contains too many nonzero entries %i %i", hmfile, line, count+1, MAXLENINDEXMODE);
    			  exit(1);
    		  }
    	  }

      }
      else {
    	  if (fields[n] != 0.0){
    		  if ( count < MAXLENINDEXMODE ){
    		  eigl_mono[currmode][count] = 3*(fields[n]-1);
    		  val_eigl_mono[currmode][count] = 1.0;
    		  eigl_mono[currmode+1][count] = 3*(fields[n]-1)+1;
    		  val_eigl_mono[currmode+1][count] = 1.0;
    		  eigl_mono[currmode+2][count] = 3*(fields[n]-1)+2;
    		  val_eigl_mono[currmode+2][count] = 1.0;
    		  count ++;
    		  }
    		  else{
    			  fprintf(stderr, "Reading error in %s, line %d index mode contains more than %i  nonzero entries", hmfile, line, MAXLENINDEXMODE);
    			  exit(1);
    		  }
    	  }
      }
    }
 /*   for (int c=0; c< MAXLENINDEXMODE; c++){
    	fprintf(stderr,"%i %i %i %f\n", currlig-1, currmode-1, eigl_multi[c][currmode-1][currlig-1],val_eigl_multi[c][currmode-1][currlig-1]);
    }*/
    currmode += 3;
  }

  fprintf(stderr,"%d index modes read for %s %d\n", currmode, hmword, currlig);
  nhm[currlig-1] = currmode;
  if (currlig != nlig) {
    fprintf(stderr, "Reading error in %s, line %d: Read %d %ss, expected %d ligands\n", hmfile, line, currlig, hmword, nlig);
    exit(1);
  }


}
