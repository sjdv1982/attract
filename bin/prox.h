#ifndef PROX_H
#define PROX_H

#include "max.h"
#include <cmath>

/*
const double proxlim = 16;
const double proxmax = 200;
const int proxspace = 2000;
const int proxconst = 800000;
const int proxarsize = ceil((1/proxlim-1/proxmax)*proxconst); //46000
*/

/*const double proxlim = 36;*/
/*const double proxmax = 200;*/
const int proxspace = 2000;
const int proxconst = 800000;


extern int *proxmap;
extern bool proxmap_initialized;

struct Prox {
  double plateaudissq;
  double proxlim;
  double proxmax;
  int proxmaxtype;
  double *prox[MAXATOMTYPES][MAXATOMTYPES]; 
};  

extern Prox *prox_init(int cartstatehandle, double plateaudissq, double proxlim, double proxmax, int proxmaxtype, bool has_pot);

#endif
