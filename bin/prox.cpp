#include "prox.h"
#include "state.h"

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>

#include <cstring>

int *proxmap = 0;
bool proxmap_initialized = 0;

CartState &cartstate_get(int handle);

void proxmapinit(double proxlim, double proxmax) {
  double proxliminv = 1.0/proxlim;
  int proxmapsize = ceil((proxmax-proxlim)*proxspace);
  proxmap = new int[proxmapsize];
  for (int n = 0; n < proxmapsize; n++) {
    double dsq = proxlim + double(n)/proxspace;
    double ddsq = 1.0/dsq;
    int pos = int((proxliminv-ddsq)*proxconst+0.5);
    proxmap[n] = pos;
  }
}

Prox *proxes[100];
int proxcount = 0;
char *shmlinks[100];
int shmlinkcount = 0;

Prox *prox_init(int cartstatehandle, double plateaudissq, double proxlim, double proxmax, int proxmaxtype, bool has_pot) {
  
  const int proxarsize = ceil((1/proxlim-1/proxmax)*proxconst);
  Prox *p;
  for (int n = 0; n < proxcount; n++) {
    p = proxes[n];
    if (fabs(p->plateaudissq - plateaudissq) < 0.01
     &&  fabs(p->proxlim - proxlim) < 0.01
     &&  fabs(p->proxmax - proxmax) < 0.01
     &&  abs(p->proxmaxtype - proxmaxtype) == 0
     ) return p;
  } 
  
  
  proxes[proxcount] = new Prox;
  
  p = proxes[proxcount];
  proxcount++;
  p->plateaudissq = plateaudissq;
  p->proxlim = proxlim;
  p->proxmax = proxmax;
  p->proxmaxtype = proxmaxtype;
  
  CartState &cartstate = cartstate_get(cartstatehandle); 

  if (!proxmap_initialized) {
    proxmapinit(proxlim, proxmax);
    proxmap_initialized = 1;
  }
  double proxliminv = 1/proxlim;
  int proxmaxcount = 0;
  for (int i = 0; i < proxmaxtype; i++) {
    for (int j = i; j < proxmaxtype; j++) {
      proxmaxcount++;
    }
  }

  char shm_proxname[100];
  sprintf(shm_proxname, "/attract-prox-%d-%d-%d-%d", int(1000*plateaudissq),int(1000*proxlim),int(1000*proxmax), proxmaxtype);     
  int fshmprox = shm_open(shm_proxname, O_RDONLY, S_IREAD);
  
  double *proxdata;
  bool loaded = 0;
  int proxdatasize = proxmaxcount*2*proxarsize*sizeof(double);
  if (fshmprox != -1) {
    ftruncate(fshmprox, proxdatasize);
    proxdata = (double *) mmap(0,proxdatasize,
     PROT_READ, (MAP_SHARED | MAP_NORESERVE), fshmprox, 0);
    if (proxdata == MAP_FAILED) {
      fprintf(stderr, "Error in loading shared prox map\n");
      exit(0);
    } 
    loaded = 1;
  }
  else {
    fshmprox = shm_open(shm_proxname, (O_CREAT| O_RDWR), (S_IREAD | S_IWRITE));
    if (fshmprox == -1) {
      fprintf(stderr, "Error in creating shared prox map\n");
      exit(0);      
    }
    char *shmlinkname = new char[100];
    strcpy(shmlinkname, shm_proxname);
    shmlinks[shmlinkcount] = shmlinkname;
    shmlinkcount++;
    ftruncate(fshmprox, proxdatasize);
    proxdata = (double *) mmap(0,proxdatasize,
     (PROT_READ | PROT_WRITE), MAP_SHARED, fshmprox, 0);
    if (proxdata == MAP_FAILED) {
      fprintf(stderr, "Error in creating shared prox map (map error)\n");
      perror("");
      shm_unlink(shmlinkname);
      exit(0);
    } 
    memset(proxdata, 0, proxdatasize);
  }
  
  int proxcount = 0;
  float swi_on = cartstate.swi_on;  
  float swi_off = cartstate.swi_off;
  int potshape =  cartstate.potshape;
  for (int i = 0; i < proxmaxtype; i++) {
    for (int j = i; j < proxmaxtype; j++) {
      //TODO: invoke nonbon.h
      double *currprox = &proxdata[proxcount];
      p->prox[i][j] = currprox;
      p->prox[j][i] = currprox;
      proxcount += 2*proxarsize;
      if (loaded) continue;

      double rci = cartstate.rc[i][j];
      double aci = cartstate.ac[i][j];
      double emini = cartstate.emin[i][j];
      double rmin2i = cartstate.rmin2[i][j];
      int ivor = cartstate.ipon[i][j];      
      for (int n = 0; n < proxarsize;n++) {
	double rr2 = proxliminv - double(n)/proxconst;
	double dsq = 1.0/rr2;
	double alen = aci;
	double rlen = rci;

	double rr23 = rr2 * rr2 * rr2;
        double rrd;
	if (potshape==8) {
	  rrd = rr2;
	}
	else if (potshape==12) {
	  rrd = rr23;
	}	
			
	double rep = rlen * rrd;
	double vlj = (rep-alen)*rr23; 
	double fb=6.0*vlj+2.0*(rep*rr23);	

	double fswi = 1;
	if (swi_on > 0 || swi_off > 0) {
	  if (dsq > swi_on*swi_on) {
	    if (dsq > swi_off*swi_off) {
	      fswi = 0;
	    }
	    else {
	      double distance = sqrt(dsq) ;
	      fswi = 1-(distance - swi_on)/(swi_off-swi_on);
	    }
	  }    
	}

	double energy, grad;
	if (dsq < rmin2i) {
	  energy = fswi * (vlj + (ivor-1) * emini);
	  grad = fswi * (fb * rr2);
	}
	else {
	  energy = fswi * ivor * vlj;
	  grad = fswi * ivor * fb * rr2;
	}
        if (has_pot) {
	  rr2 = 1.0/plateaudissq;
	  rr23 = rr2 * rr2 * rr2;
	  if (potshape==8) {
	    rrd = rr2;
	  }
	  else if (potshape==12) {
	    rrd = rr23;
	  }

	  rep = rlen * rrd;
	  vlj = (rep-alen)*rr23; 
	  fb=6.0*vlj+2.0*(rep*rr23);	

          fswi = 1;
	  if (swi_on > 0 || swi_off > 0) {
	    if (plateaudissq > swi_on*swi_on) {
	      if (plateaudissq > swi_off*swi_off) {
		fswi = 0;
	      }
	      else {
		double distance = sqrt(plateaudissq) ;
		fswi = 1-(distance - swi_on)/(swi_off-swi_on);
	      }
	    }    
	  }

	  if (plateaudissq < rmin2i) {
	    energy -= fswi * (vlj + (ivor-1) * emini);
	    grad -= fswi * fb * rr2 * sqrt(dsq/plateaudissq);
	  }
	  else {
	    energy -= fswi * ivor * vlj;
	    grad -= fswi * ivor * fb * rr2 * sqrt(dsq/plateaudissq);
	  }
	}
	currprox[2*n] = energy;
	currprox[2*n+1] = grad;	
      }
    }
  }
  return p;
}

