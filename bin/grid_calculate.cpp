#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>

#include "state.h"
#include "makegrid.h"
#include "nonbon.h"

extern "C" void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi);

extern CartState &cartstate_get(int handle);

inline void inc(Coor *dis, int nrdis, int index, double d) {
  for (int n = 0; n < nrdis; n++) 
    dis[n][index] += d;
};

int ener_max = 0;
inline unsigned int new_energrad(Grid &g) {
  if (g.nr_energrads == ener_max) {
    int new_ener_max = 2 * ener_max;
    if (ener_max == 0) new_ener_max = 100000;
    g.energrads = (EnerGrad *) realloc(g.energrads,
      new_ener_max*sizeof(EnerGrad));
    memset(g.energrads+ener_max,0,(new_ener_max-ener_max)*sizeof(EnerGrad));
    ener_max = new_ener_max;
  }
  g.nr_energrads++;
  return g.nr_energrads;
}

void Grid::calculate(int cartstatehandle, int ligand, const char *interior_grid, double plateaudis, double neighbourdis, int gridextension, int nhm0, bool (&alphabet)[MAXATOMTYPES]) {  
    
  init(gridspacing, gridextension, plateaudis, neighbourdis, alphabet);
  neighbours = new Neighbour[100000000]; //max 100 million neighbours
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
  
  //Initialize grids
  innergrid = new Voxel[gridx*gridy*gridz];
  memset(innergrid, 0, gridx*gridy*gridz*sizeof(Voxel));

  gridx2 = int((gridx + 2 * gridextension)/2)+1;
  gridy2 = int((gridy + 2 * gridextension)/2)+1;
  gridz2 = int((gridz + 2 * gridextension)/2)+1;
  biggrid = new Potential[gridx2*gridy2*gridz2];
  memset(biggrid, 0, gridx2*gridy2*gridz2*sizeof(Potential));
  
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
  
  Coor *disx = dis; //grid distances
  inc(disx, nratoms,  0, -gridextension*gridspacing);
  inc(disx, nratoms,  1, -gridextension*gridspacing);
  inc(disx, nratoms,  2, -gridextension*gridspacing);

  Coor *disjx = new Coor[nratoms]; 
  memcpy(disjx, disx, nratoms * sizeof(Coor));
  inc(disjx, nratoms,  0, 0.5*gridspacing);
  inc(disjx, nratoms,  1, 0.5*gridspacing);
  inc(disjx, nratoms,  2, 0.5*gridspacing);
 
  Coor *disxy = new Coor[nratoms];
  Coor *disxyz = new Coor[nratoms];

  Coor *disjxy = new Coor[nratoms];
  Coor *disjxyz = new Coor[nratoms];
     
  bool inner_x, inner_xy, inner_xyz;  
  bool junction_x, junction_xy, junction_xyz;
  int outerminx = 2*(int((gridx - 1)/2) - 1);
  int outerminy = 2*(int((gridy - 1)/2) - 1);
  int outerminz = 2*(int((gridz - 1)/2) - 1);
  
  
  //Main loop
  int potcount = 0;
  for (int x = -gridextension; x < gridx+gridextension; x++) {    
    fprintf(stderr, "%d/%d %d %d\n", x, gridx+gridextension, potcount, nr_neighbours);  
        
    memcpy(disxy, disx, nratoms*sizeof(Coor));
    memcpy(disjxy, disjx, nratoms*sizeof(Coor));
    inner_x = 1;    
    if (x < 0 || x >= gridx) inner_x = 0;
    junction_x = 1;
    if (x < 0 || x >= gridx-1) junction_x = 0;
        
    for (int y = -gridextension; y < gridy+gridextension; y++) {
      inner_xy = inner_x;    
      if (y < 0 || y >= gridy) inner_xy = 0;
      junction_xy = junction_x;
      if (y < 0 || y >= gridy-1) junction_xy = 0;
        
      memcpy(disxyz, disxy, nratoms*sizeof(Coor));    
      memcpy(disjxyz, disjxy, nratoms*sizeof(Coor));    
      
      for (int z = -gridextension; z < gridz+gridextension; z++) {
	 inner_xyz = inner_xy;    
	 if (z < 0 || z >= gridz) inner_xyz = 0;
	 junction_xyz = junction_xy;
	 if (z < 0 || z >= gridz-1) junction_xyz = 0;
	 bool inside = 0;
	 if (inner_xyz) { 
	   long index = x + gridx * y + gridx*gridy*z;
	   inside = (fabs(grid[index]-interior_value) < 0.1);
   	   if (!inside) { 
	     Voxel &v = innergrid[index];
	     potcount++;
	     _calc_potential(v.potential, disxyz, xb, nratoms, wer,charges,atomtypes,
	      rc,ac,emin,rmin2,ipon,potshape); 
	     if (junction_xyz) {
	       _calc_neighbours(
	        v.neighbourlist, v.nr_neighbours, disjxyz, nratoms,atomtypes
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
	   int zz = (z + gridextension)/2;
	   if (xx < 0 || xx >= gridx2) printf("ERR %d %d\n", xx, gridx2);
	   if (yy < 0 || yy >= gridy2) printf("ERR %d %d\n", yy, gridy2);
	   if (zz < 0 || zz >= gridz2) printf("ERR %d %d\n", zz, gridz2);

	   long index = xx + gridx2 * yy + gridx2*gridy2*zz;	   
	   Potential &p = biggrid[index];
	   if (p[MAXATOMTYPES]) printf("ERR %d %d %d\n", xx, yy, zz);
           potcount++;	   
           _calc_potential(p, disxyz, xb, nratoms, wer,charges,atomtypes,
	    rc,ac,emin,rmin2,ipon,potshape); 	   
	 }
	 inc(disxyz, nratoms,  2, gridspacing);
	 inc(disjxyz, nratoms,  2, gridspacing);
      }
      
      inc(disxy, nratoms,1, gridspacing);
      inc(disjxy, nratoms,1, gridspacing);
    }
      
    inc(disx, nratoms,0, gridspacing);
    inc(disjx, nratoms,0, gridspacing);
  }     
    
}
#ifdef TORQUEGRID
inline void Grid::_calc_potential_elec(float &energy, Coorf &grad, float (&torques)[9], const  Coor *dis, const Coor *xb, int nrdis, double *charges) {
#else
inline void Grid::_calc_potential_elec(float &energy, Coorf &grad, 
const  Coor *dis, const Coor *xb, int nrdis, double *charges) {
#endif
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

inline void Grid::_calc_potential(Potential &p, const Coor *dis, const Coor *xb, int nrdis,
double *wer, double *charges,int *atomtypes,
const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, const int &potshape //ATTRACT params
) {
  memset(p, 0, sizeof(Potential));
  for (int i = 0; i <= MAXATOMTYPES; i++) {
    if (i == MAXATOMTYPES || alphabet[i]) {
      p[i] = new_energrad(*this);
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
      
      nonbon(1,wer[n],rci,aci,emini,rmin2i,ivor,dsq,
       1/dsq, dd[0],dd[1],dd[2],potshape,
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
#ifdef TORQUEGRID  
  _calc_potential_elec(e.energy, e.grad, 
  e.torques, dis, xb, nrdis, charges);  
#else  
  _calc_potential_elec(e.energy, e.grad, 
  dis, xb, nrdis, charges);  
#endif
}

inline void Grid::_calc_neighbours(int &neighbourlist, short &neighboursize, const Coor *dis, int nrdis, int *atomtypes) {
  for (int n = 0; n < nrdis; n++) {
    if (atomtypes[n] == 0) continue;  
    const Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];       
    if (dsq >= neighbourdissq) continue;
    Neighbour &nb = neighbours[nr_neighbours];
    nr_neighbours++;
    if (neighboursize == 0) {
      neighbourlist = nr_neighbours;
    }
    neighboursize++;
    nb.type = 2;
    if (dsq <= plateaudissq) nb.type = 1;    
    nb.index = n;
  }  
}
