#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>

#include "state.h"
#include "makegrid.h"
#include "nonbon.h"

int get_new_task(int *taskstat, int nrtasks) {
  //returns task, -1 if no free tasks, -2 if done;
  for (int n = 0; n < nrtasks; n++) {
    if (taskstat[n] == 0) {
      taskstat[n] = 1;
      return n;
    }
  }
  for (int n = 0; n < nrtasks; n++) {
    if (taskstat[n] == 1) return -1;
  }
  return -2;  
}

inline void _calc_potential_elec(float &energy, Coorf &grad, 
#ifdef TORQUEGRID
float (&torques)[9],
#endif
const  Coor *dis, const Coor *xb, int nrdis, double *charges,
double plateaudissq
);

inline void _calc_potential(Potential &p, const Coor *dis, const Coor *xb, 
 int nrdis, double *wer, double *charges, int *atomtypes,
 const Parameters &rc, const Parameters &ac, 
 const Parameters &emin, const Parameters &rmin2, const iParameters &ipon, 
 const int &potshape, const float &swi_on, const float &swi_off,
 bool (&alphabet)[MAXATOMTYPES],
 int &nr_energrads, EnerGrad *&energrads,
 double plateaudissq
);
inline short _calc_neighbours(int &neighbourlist, 
 const Coor *dis, int nrdis, int *atomtypes,
 int neighbourdissq, double plateaudissq,
 int &nr_neighbours, Neighbour *neighbours, Neighbour *neigh
);

extern "C" void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi);

extern CartState &cartstate_get(int handle);

inline void inc(Coor *dis, int nrdis, int index, double d) {
  for (int n = 0; n < nrdis; n++) 
    dis[n][index] += d;
};

int ener_max = 0;
inline int new_energrad(int &nr_energrads, EnerGrad *&energrads) {  
  int ret;
#pragma omp critical
{    
  if (nr_energrads == ener_max) {
    int new_ener_max = 2 * ener_max;
    if (ener_max == 0) new_ener_max = 100000;
    energrads = (EnerGrad *) realloc(energrads,
      new_ener_max*sizeof(EnerGrad));
    memset(energrads+ener_max,0,(new_ener_max-ener_max)*sizeof(EnerGrad));
    ener_max = new_ener_max;
  }
  nr_energrads++;
  ret = nr_energrads;
}
  return ret;
}

void calculate_line(
  int gridextension,
  int gridx,
  int gridy,
  int gridz,
  int nratoms,
  const Coor *xb,
  double neighbourdissq, double plateaudissq,
  bool calc_pot,
  
  int z,
  Coor *dis,
  
  double *wer,  double *charges, int *atomtypes,
  const Parameters &rc, const Parameters &ac, 
  const Parameters &emin, const Parameters &rmin2, const iParameters &ipon,   
  int potshape, float swi_on, float swi_off,
  bool (&alphabet)[MAXATOMTYPES],
  
  double *cgrid,
  Voxel *cinnergrid,
  Potential *cbiggrid,
  int &potcount,
  Neighbour *neighbours,
  int &nr_neighbours,
  EnerGrad *&energrads,
  int &nr_energrads
) 
{  
  Coor *disz = new Coor[nratoms]; 
  memcpy(disz, dis, nratoms * sizeof(Coor));
  inc(disz, nratoms,  0, -gridextension*gridspacing);
  inc(disz, nratoms,  1, -gridextension*gridspacing);
  inc(disz, nratoms,  2, z*gridspacing);

  Coor *disjz = new Coor[nratoms]; 
  memcpy(disjz, disz, nratoms * sizeof(Coor));
  inc(disjz, nratoms,  0, 0.5*gridspacing);
  inc(disjz, nratoms,  1, 0.5*gridspacing);
  inc(disjz, nratoms,  2, 0.5*gridspacing);

  bool inner_z = 1;    
  if (z < 0 || z >= gridz) inner_z = 0;
  bool junction_z = 1;
  if (z < 0 || z >= gridz-1) junction_z = 0;

  Coor *diszy = new Coor[nratoms];
  Coor *diszyx = new Coor[nratoms];

  Coor *disjzy = new Coor[nratoms];
  Coor *disjzyx = new Coor[nratoms];

  memcpy(diszy, disz, nratoms*sizeof(Coor));
  memcpy(disjzy, disjz, nratoms*sizeof(Coor));

  int gridx2 = int((gridx + 2 * gridextension)/2)+1;
  int outerminx = 2*(int((gridx - 1)/2) - 1);
  int outerminy = 2*(int((gridy - 1)/2) - 1);
  int outerminz = 2*(int((gridz - 1)/2) - 1);

  bool inner_zy, inner_zyx;
  bool junction_zy, junction_zyx;

  Neighbour *neigh = new Neighbour[32768];  
  memset(neigh, 0, 32768*sizeof(Neighbour));
  
  for (int y = -gridextension; y < gridy+gridextension; y++) {
    inner_zy = inner_z;    
    if (y < 0 || y >= gridy) inner_zy = 0;
    junction_zy = junction_z;
    if (y < 0 || y >= gridy-1) junction_zy = 0;

    memcpy(diszyx, diszy, nratoms*sizeof(Coor));    
    memcpy(disjzyx, disjzy, nratoms*sizeof(Coor));    
    for (int x = -gridextension; x < gridx+gridextension; x++) {
       inner_zyx = inner_zy;    
       if (x < 0 || x >= gridx) inner_zyx = 0;
       junction_zyx = junction_zy;
       if (x < 0 || x >= gridx-1) junction_zyx = 0;
       bool inside = 0;
       if (inner_zyx) { 
	 int index = x + gridx * y;
	 inside = (fabs(cgrid[index]-interior_value) < 0.1);
   	 if (!inside) { 
	   Voxel &v = cinnergrid[index];
	   if (calc_pot) {
	     #pragma omp atomic
	     potcount++;
	     _calc_potential(v.potential, diszyx, xb, nratoms, 
	      wer,charges,atomtypes,
	      rc,ac,emin,rmin2,ipon,potshape,swi_on, swi_off,
	      alphabet,nr_energrads,energrads,plateaudissq
	     ); 
	   }
	   if (junction_zyx) {
	     v.nr_neighbours = _calc_neighbours(
	      v.neighbourlist, disjzyx, nratoms,atomtypes,
	      neighbourdissq,plateaudissq,nr_neighbours, neighbours, neigh
	     );
	   } 
	 }	   
       }
       if (!(x % 2) && !(y % 2) && !(z % 2) &&
       (
       ((x <= 0) || (x >= outerminx)) ||
       ((y <= 0) || (y >= outerminy)) ||
       ((z <= 0) || (z >= outerminz)) 
       )
       )
       {
	 int xx = (x + gridextension)/2;
	 int yy = (y + gridextension)/2;	   

	 int index = xx + gridx2 * yy;
	 if (calc_pot) {
	   Potential &p = cbiggrid[index];
	   #pragma omp atomic
           potcount++;	   
           _calc_potential(p, diszyx, xb, nratoms, wer,charges,atomtypes,
	    rc,ac,emin,rmin2,ipon,potshape, swi_on, swi_off,
	    alphabet,nr_energrads,energrads,plateaudissq
	   ); 	   
	 }
       }
       
       inc(diszyx, nratoms,  0, gridspacing);
       inc(disjzyx, nratoms,  0, gridspacing);
    }
    inc(diszy, nratoms,1, gridspacing);
    inc(disjzy, nratoms,1, gridspacing);
  }

  delete[] diszy;
  delete[] diszyx;
  delete[] disjzy;
  delete[] disjzyx;
  
  delete[] disz;
  delete[] disjz;
  delete[] neigh;
}

#define MACRO_calculate_line \
    calculate_line( \
     gridextension, \
     gridx, \
     gridy, \
     gridz, \
     nratoms, \
     xb, \
     neighbourdissq, plateaudissq, \
     calc_pot, \
     \
     z, \
     dis, \
     \
     wer,  charges, atomtypes, \
     rc, ac,  \
     emin, rmin2, ipon, \
     potshape, swi_on, swi_off, \
     alphabet, \
     \
     cgrid, \
     cinnergrid, \
     cbiggrid, \
     potcount, \
     neighbours, \
     nr_neighbours, \
     energrads, \
     nr_energrads \
    ); 


void Grid::calculate(int cartstatehandle, int ligand, const char *interior_grid, double plateaudis, double neighbourdis, int gridextension, int nhm0, bool (&alphabet)[MAXATOMTYPES], bool calc_pot) {  
  init(gridspacing, gridextension, plateaudis, neighbourdis, alphabet);

  nr_neighbours = 0; 
  shm_neighbours = -1; 
  energrads = NULL;
  nr_energrads = 0;
  shm_energrads = -1;
  nhm = nhm0;  
  
  int n;
  char interior_grid2[1000];
  if (strlen(interior_grid) > 1000) {
    fprintf(stderr, "Name of grid too large\n");
    exit(1);
  }
  strcpy(interior_grid2, interior_grid);
  
  //Load interior grid
  unsigned int gridx0, gridy0, gridz0;
  double *grid;
  double gridspacing0;
  read_vol(interior_grid2, &gridspacing0, &ori[0], &ori[1],&ori[2],&gridx0,&gridy0,&gridz0,&grid);
  
  gridspacing = gridspacing0;
  gridx = gridx0; gridy = gridy0; gridz = gridz0;
  
  //Load Cartesian state and parameters
  CartState &cartstate = cartstate_get(cartstatehandle); 

  int nratoms = cartstate.natom[ligand];
  int start = 0;
  if (ligand > 0) start = cartstate.ieins[ligand-1];
  double *x = &(cartstate.x[3*start]);
  pivot[0] = cartstate.pivot[0][ligand];
  pivot[1] = cartstate.pivot[1][ligand];
  pivot[2] = cartstate.pivot[2][ligand];
  Coor *xb = (Coor *) &(cartstate.xb[0]);
  double *charges = &(cartstate.chai[start]);
  int *atomtypes = &(cartstate.iaci[start]);
  double *wer = &(cartstate.we[start]);
  Parameters &rc = cartstate.rc;
  Parameters &ac = cartstate.ac;
  Parameters &emin = cartstate.emin;
  Parameters &rmin2 = cartstate.rmin2;
  iParameters &ipon = cartstate.ipon;
  int potshape = cartstate.potshape;
  float swi_on = cartstate.swi_on;
  float swi_off = cartstate.swi_off;
    
  //Initialize grids
  gridx2 = int((gridx + 2 * gridextension)/2)+1;
  gridy2 = int((gridy + 2 * gridextension)/2)+1;
  gridz2 = int((gridz + 2 * gridextension)/2)+1;

  Voxel *cinnergrid;
  Potential *cbiggrid;
  
  innergrid = new Voxel[gridx*gridy*gridz];
  memset(innergrid, 0, gridx*gridy*gridz*sizeof(Voxel));

  if (calc_pot) {
    biggrid = new Potential[gridx2*gridy2*gridz2];
    memset(biggrid, 0, gridx2*gridy2*gridz2*sizeof(Potential));
  }
  neighbours = new Neighbour[100000000]; //max 100 million neighbours
  
  //Set up distances 
  Coor *dis0 = new Coor[nratoms]; 
  for (n = 0; n < nratoms; n++) {
    dis0[n][0] = -(x[3*n]-ori[0]);
    dis0[n][1] = -(x[3*n+1]-ori[1]);
    dis0[n][2] = -(x[3*n+2]-ori[2]);
  }
  Coor *dis = new Coor[nratoms]; 
  natoms = nratoms;
  memcpy(dis, dis0, nratoms * sizeof(Coor));
          
  //Main loop
  int potcount = 0;
  #pragma omp parallel for private(cinnergrid,cbiggrid) schedule(dynamic,1) 
  for (int z = -gridextension; z < gridz+gridextension; z++) {    
    fprintf(stderr, "%d/%d %d %d\n", z, gridz+gridextension, potcount, nr_neighbours);  
        
    double *cgrid = grid + gridx*gridy*z;
    cinnergrid = innergrid + gridx*gridy*z;
    int zz = (z + gridextension)/2;
    cbiggrid = biggrid + gridx2*gridy2*zz;
    
    MACRO_calculate_line
  }     
  fprintf(stderr, "%d/%d %d %d\n", gridz+gridextension, gridz+gridextension, potcount, nr_neighbours);    

}
inline void _calc_potential_elec(float &energy, Coorf &grad, 
#ifdef TORQUEGRID
float (&torques)[9],
#endif
const  Coor *dis, const Coor *xb, int nrdis, double *charges,
double plateaudissq

) {
  double energy0 = 0; Coor grad0 = {0,0,0};
  for (int n = 0; n < nrdis; n++) {
    float charge = charges[n];
    if (fabs(charge) > 0.001) {    
      const Coor &d = dis[n];
      double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
      double rr2 = 1.0/dsq;
      double ratio = 1;
      if (dsq < plateaudissq) {
	ratio = sqrt(dsq/plateaudissq);
	dsq = plateaudissq;
	rr2 = 1.0/dsq;
      }
      Coor dd = {d[0]*rr2,d[1]*rr2,d[2]*rr2};          
      
      double charge0 = charge * felecsqrt; //at realtime, to be multiplied with the charge of the other atom...
      
      double energy00; Coor grad00;
      elec(1,charge0,rr2,dd[0]*ratio,dd[1]*ratio,dd[2]*ratio,energy00,grad00);
      energy0 += energy00;
      grad0[0] += grad00[0];
      grad0[1] += grad00[1];
      grad0[2] += grad00[2];
#ifdef TORQUEGRID       
      torques[0] -= xb[n][0] * grad00[0];
      torques[1] -= xb[n][0] * grad00[1];
      torques[2] -= xb[n][0] * grad00[2];
      torques[3] -= xb[n][1] * grad00[0];
      torques[4] -= xb[n][1] * grad00[1];
      torques[5] -= xb[n][1] * grad00[2];
      torques[6] -= xb[n][2] * grad00[0];
      torques[7] -= xb[n][2] * grad00[1];
      torques[8] -= xb[n][2] * grad00[2];
#endif      
    }
  }
  energy = energy0;
  grad[0] = grad0[0];
  grad[1] = grad0[1];
  grad[2] = grad0[2];
}

inline void _calc_potential(Potential &p, const Coor *dis, const Coor *xb, int nrdis,
double *wer, double *charges,int *atomtypes,
const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, const int &potshape, const float &swi_on, const float &swi_off, //ATTRACT params

bool (&alphabet)[MAXATOMTYPES],
int &nr_energrads, EnerGrad *&energrads,
double plateaudissq //grid params
) {
  memset(&p, 0, sizeof(Potential));
  for (int i = 0; i <= MAXATOMTYPES; i++) {
    if (i == MAXATOMTYPES || alphabet[i]) {
      p[i] = new_energrad(nr_energrads, energrads);
    }        
  }
  for (int n = 0; n < nrdis; n++) {
    int atomtype = atomtypes[n];
    if (atomtype == 0) continue;  
    const Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2]; 
    if (dsq >= distcutoffsq) continue;

    double rr2 = 1.0/dsq;
    double ratio = 1;
    if (dsq < plateaudissq) {
      ratio = sqrt(dsq/plateaudissq);
      dsq = plateaudissq;
      rr2 = 1.0/dsq;
    }
    Coor dd = {d[0]*rr2*ratio,d[1]*rr2*ratio,d[2]*rr2*ratio};          

    for (int i = 0; i < MAXATOMTYPES; i++) {
      if (!alphabet[i]) continue;
      double rci = rc[atomtype-1][i];
      double aci = ac[atomtype-1][i];
      double emini = emin[atomtype-1][i];
      double rmin2i = rmin2[atomtype-1][i];
      int ivor = ipon[atomtype-1][i];
      
      double energy; Coor grad;
      
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
      
      nonbon(1,wer[n],rci,aci,emini,rmin2i,ivor,dsq,
       1/dsq, dd[0],dd[1],dd[2],potshape, fswi,
       energy, grad);
      EnerGrad &e = *(energrads + p[i] - 1);
      e.energy += energy;
      e.grad[0] += grad[0];
      e.grad[1] += grad[1];
      e.grad[2] += grad[2];
#ifdef TORQUEGRID      
      e.torques[0] -= xb[n][0] * grad[0];
      e.torques[1] -= xb[n][0] * grad[1];
      e.torques[2] -= xb[n][0] * grad[2];
      e.torques[3] -= xb[n][1] * grad[0];
      e.torques[4] -= xb[n][1] * grad[1];
      e.torques[5] -= xb[n][1] * grad[2];
      e.torques[6] -= xb[n][2] * grad[0];
      e.torques[7] -= xb[n][2] * grad[1];
      e.torques[8] -= xb[n][2] * grad[2];
#endif      
    }
  }
  EnerGrad &e = *(energrads + p[MAXATOMTYPES] - 1);
  _calc_potential_elec(e.energy, e.grad, 
#ifdef TORQUEGRID    
  e.torques,
#endif
  dis, xb, nrdis, charges, plateaudissq);  
}

inline short _calc_neighbours(int &neighbourlist, const Coor *dis, int nrdis, int *atomtypes, 
int neighbourdissq, double plateaudissq,
int &nr_neighbours, Neighbour *neighbours, Neighbour *neigh
) {
  short neighboursize = 0;
  for (int n = 0; n < nrdis; n++) {
    if (atomtypes[n] == 0) continue;  
    const Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];       
    if (dsq >= neighbourdissq) continue;
    Neighbour &nb = neigh[neighboursize];
    neighboursize++;
    if (neighboursize == 32766) fprintf(stderr, "Neighbour overflow\n");
    nb.type = 2;
    if (dsq <= plateaudissq) nb.type = 1;    
    nb.index = n;
  }

  if (neighboursize) {
    #pragma omp critical
{  
    neighbourlist = nr_neighbours+1;
    memcpy(neighbours+nr_neighbours, neigh, neighboursize*sizeof(Neighbour));
    nr_neighbours += neighboursize;
}
  }
  
  return neighboursize;
}
