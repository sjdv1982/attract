#include "grid.h"
#include "nonbon.h"
#include "state.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

void Grid::init(double gridspacing0, int gridextension0, double plateaudis0,
double neighbourdis0, bool (&alphabet0)[MAXATOMTYPES]) {
  gridspacing = gridspacing0;
  gridextension = gridextension0;
  plateaudis = plateaudis0;
  plateaudissq = plateaudis*plateaudis;
  plateaudissqinv = 1.0/plateaudissq;
  neighbourdis = neighbourdis0;
  neighbourdissq = neighbourdis * neighbourdis;
  architecture = ARCHITECTURE;
  //Pre-compute the scale-down-distance ratios
  int size_ratio  = int(10000*plateaudissq);
  _ratio = new double[size_ratio+1];
  for (int n = 0; n <= size_ratio; n++) {
    double dissq = ((n+0.5)/10000);
    _ratio[n] = sqrt(dissq/plateaudissq);
  }
  memset(alphabet, 0, MAXATOMTYPES*sizeof(bool));
  memcpy(alphabet, alphabet0, sizeof(alphabet));
  alphabetsize = 0;
  for (int n = 0; n < MAXATOMTYPES; n++) {
    if (alphabet[n]) {
      alphabetsize++;
    }
  }
} 

extern "C" void nonbon_grid_std(
  const Grid *&g, const int &rigid, 
  const int &iab, const int &fixre,
  const Coor *xl, const Coor *xr,const Coor &pivotr,const Coor &tr,  
  const double *wel, const double *wer, const double *chail, const double *chair, const int *iacil, const int *iacir, 
  const int &natoml,const int &natomr,

const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, const int &potshape, const int &cdie, const double &epsilon,
  const float &swi_on, const float &swi_off, 
  //ATTRACT params
  
  Coor *fl, double &evdw, double &eelec,
  
  Coor *fr, const double (&pm2)[3][3][3], double *deltar);



extern "C" void nonbon_grid_torque(
  const Grid *&g, const int &rigid, 
  const int &iab, const int &fixre,
  const Coor *xl, const Coor *xr,const Coor &pivotr,const Coor &tr,  
  const double *wel, const double *wer, const double *chail, const double *chair, const int *iacil, const int *iacir, 
  const int &natoml,const int &natomr,

const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, const int &potshape, const int &cdie, const double &epsilon,
  const float &swi_on, const float &swi_off, 
  //ATTRACT params
  
  Coor *fl, double &evdw, double &eelec,
  
  Coor *fr, const double (&pm2)[3][3][3], double *deltar);

void get_shm_name(int shm_id, char *shm_name) {
  sprintf(shm_name, "/attract-grid%d", shm_id);
}

extern "C" void nonbon_grid_(
  const Grid *&g, const int &torquegrid, const int &rigid, 
  const int &iab, const int &fixre,
  const Coor *xl, const Coor *xr,const Coor &pivotr,const Coor &tr,  
  const double *wel, const double *wer, const double *chail, const double *chair, const int *iacil, const int *iacir, 
  const int &natoml,const int &natomr,

const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, const int &potshape, const int &cdie, const double &epsilon,
  const float &swi_on, const float &swi_off, 
  //ATTRACT params
  
  Coor *fl, double &evdw, double &eelec,
  
  Coor *fr, const double (&pm2)[3][3][3], double *deltar)
{
  if (torquegrid) {
    nonbon_grid_torque(
     g,rigid,iab,fixre,xl,xr,pivotr,tr,wel,wer,chail,chair,iacil,iacir,natoml,natomr,
     rc,ac,emin,rmin2,ipon,potshape,cdie,epsilon,swi_on,swi_off,
     fl,evdw,eelec,fr,pm2,deltar
    ); 
  }    
  else {  
    nonbon_grid_std(
     g,rigid,iab,fixre,xl,xr,pivotr,tr,wel,wer,chail,chair,iacil,iacir,natoml,natomr,
     rc,ac,emin,rmin2,ipon,potshape,cdie,epsilon,swi_on,swi_off,
     fl,evdw,eelec,fr,pm2,deltar
    ); 
  }
}