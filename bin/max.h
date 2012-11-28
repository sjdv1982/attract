#ifndef MAX_H 
#define MAX_H /* to make sure that we include only once... */

const int MAXSTRUC = 100000; //only for deredundant
const int MAXATOM = 5000;
const int MAXRES = 1000;
const int TOTMAXATOM = 30000;
const int TOTMAXRES = 5000;
const int MAXLIG = 7;
const int MAXMODE = 5;
/*const int MAXMOLPAIR = 25000000; //can be 5x smaller if only two-body coarse-grained docking is done*/
const int MAXMOLPAIR = 1000000; 
const int MAXDOF = 43; 
const int MAXATOMTYPES = 99;
const int MAXSELECTION = 1; //maximum size of selection; NOTE: a static array of MAXSELECTION*MAXSELECTION Coors+doubles is kept in memory!
const int MAXRESTRAINTS = 1;
const int MAXENS = 10; //maximum ensemble size

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
