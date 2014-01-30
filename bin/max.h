#ifndef MAX_H 
#define MAX_H /* to make sure that we include only once... */

const int MAXSTRUC = 200000; //only for deredundant
const int MAXATOM = 10000;
const int MAXRES = 3000;
const int TOTMAXATOM = 100000;
const int TOTMAXRES = 10000;
const int MAXLIG = 100;
const int MAXMODE = 10;
const int MAXMOLPAIR = 20000000; 
const int MAXDOF = 600; //100 ligands withouh modes
const int MAXATOMTYPES = 99;
const int MAXSELECTION = 1000; //maximum size of selection; NOTE: a static array of MAXSELECTION*MAXSELECTION Coors+doubles is kept in memory!
const int MAXRESTRAINTS = 10000;
const int MAXENS = 100; //maximum ensemble size

typedef double dof[MAXSTRUC][MAXLIG]; //only for deredundant
typedef double modes[MAXSTRUC][MAXLIG][MAXMODE]; //only for deredundant
typedef double coors[MAXSTRUC][MAXLIG][3]; //only for deredundant

typedef char char4[4];
typedef double dof2[MAXLIG];
typedef int idof2[MAXLIG];
typedef double coors2[MAXLIG][3];
typedef double modes2[MAXLIG][MAXMODE];

typedef double Parameters[MAXATOMTYPES][MAXATOMTYPES];
typedef int iParameters[MAXATOMTYPES][MAXATOMTYPES];

#endif
