#ifndef MAX_H 
#define MAX_H /* to make sure that we include only once... */

const int MAXSTRUC = 100000;
const int MAXATOM = 10000;
const int MAXRES = 500;
const int TOTMAXATOM = 100000;
const int TOTMAXRES = 5000;
const int MAXLIG = 100;
const int MAXMODE = 20;
const int MAXMOLPAIR = 5000000;
const int MAXDOF = 600; //100 ligands without modes
const int MAXATOMTYPES = 99;
const int MAXSELECTION = 1000; //maximum size of selection; NOTE: a static array of MAXSELECTION*MAXSELECTION Coors+doubles is kept in memory!
const int MAXRESTRAINTS = 10000;

typedef double dof[MAXSTRUC][MAXLIG];
typedef double modes[MAXSTRUC][MAXLIG][MAXMODE];
typedef char char4[4];
typedef double dof2[MAXLIG];
typedef double modes2[MAXLIG][MAXMODE];

typedef double Parameters[MAXATOMTYPES][MAXATOMTYPES];
typedef int iParameters[MAXATOMTYPES][MAXATOMTYPES];

#endif
