#ifndef MAX_H 
#define MAX_H /* to make sure that we include only once... */

const int MAXSTRUC = 100000;
const int MAXATOM = 30000;
const int MAXRES = 3000;
const int TOTMAXATOM = 100000;
const int TOTMAXRES = 10000;
const int MAXLIG = 50;
const int MAXMODE = 20;
const int MAXMOLPAIR = 25000000; //can be 5x smaller if only two-body coarse-grained docking is done
const int MAXDOF = 300; //50 ligands without modes
const int MAXATOMTYPES = 99;
const int MAXSELECTION = 1000; //maximum size of selection; NOTE: a static array of MAXSELECTION*MAXSELECTION Coors+doubles is kept in memory!
const int MAXRESTRAINTS = 50000;
const int MAXENS = 100; //maximum ensemble size

typedef double dof[MAXSTRUC][MAXLIG];
typedef double modes[MAXSTRUC][MAXLIG][MAXMODE];
typedef double coors[MAXSTRUC][MAXLIG][3];

typedef char char4[4];
typedef double dof2[MAXLIG];
typedef int idof2[MAXLIG];
typedef double coors2[MAXLIG][3];
typedef double modes2[MAXLIG][MAXMODE];

typedef double Parameters[MAXATOMTYPES][MAXATOMTYPES];
typedef int iParameters[MAXATOMTYPES][MAXATOMTYPES];

#endif
