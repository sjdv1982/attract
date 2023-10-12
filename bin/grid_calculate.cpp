#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>

#include "state.h"
#include "nonbon.h"

#ifdef TORQUEGRID
  #define EnerGradX EnerGradTorque
  #define energradsX energrads_torque
  #define IS_TORQUEGRID 1
#else
  #define EnerGradX EnerGradStd
  #define energradsX energrads_std
  #define IS_TORQUEGRID 0
#endif  

const double distcutoff_vdw = 50; /* cutoff for vdW interactions */
const double distcutoffsq_vdw =distcutoff_vdw * distcutoff_vdw;
const double distcutoff_elec = 50; /* cutoff for electrostatic interactions */
const double distcutoffsq_elec =distcutoff_elec * distcutoff_elec;
const int NEIGHBOUR_MAX = 32768;

static int get_new_task(int *taskstat, int nrtasks) {
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
float swi_on, float swi_off, bool cdie, double ffelec, double plateaudissq
);

inline void _calc_potential(Potential &p, const Coor *dis, const Coor *xb, 
 int nrdis, double *wer, double *charges, int *atomtypes,
 const Parameters &rc, const Parameters &ac, 
 const Parameters &emin, const Parameters &rmin2, const iParameters &ipon, 
 const int &potshape, const float &swi_on, const float &swi_off,
 int alphabetsize, bool (&alphabet)[MAXATOMTYPES],
 int &nr_energrads, EnerGradX *&energrads,
 bool cdie, double ffelec, double plateaudissq
);
inline short _calc_neighbours(int &neighbourlist, 
 const Coor *dis, int nrdis, int *atomtypes,
 double neighbourdissq, double plateaudissq,
 int &nr_neighbours, Neighbour *neighbours, Neighbour *neigh
);

extern CartState &cartstate_get(int handle);

inline void inc(Coor *dis, int nrdis, int index, double d) {
  for (int n = 0; n < nrdis; n++) 
    dis[n][index] += d;
};

static int ener_max = 0;
//(p, alphabet, energrads0, nr_energrads, energrads);
inline void new_energrads(
 Potential &p, int alphabetsize, bool (&alphabet)[MAXATOMTYPES], 
 EnerGradX (&energrads0)[MAXATOMTYPES+1], 
 int &nr_energrads, EnerGradX *&energrads  
) {
  
  EnerGradX energrads00[MAXATOMTYPES+1];
  int count = 0;
  for (int i = 0; i <= MAXATOMTYPES; i++) {
    if (i == MAXATOMTYPES || alphabet[i]) {
      memcpy(energrads00+count, energrads0+i, sizeof(EnerGradX));
      count++;
    }        
  }
  
  #pragma omp critical(energrad)
  {        
    if (nr_energrads + alphabetsize + 1 > ener_max) {
      int new_ener_max = 1.2 * ener_max;
      if (ener_max == 0) new_ener_max = 100000;
      energrads = (EnerGradX *) realloc(energrads,
        new_ener_max*sizeof(EnerGradX)
      );
      memset(energrads+ener_max,0,(new_ener_max-ener_max)*sizeof(EnerGradX));
      ener_max = new_ener_max;
    }
    memcpy(energrads + nr_energrads, energrads00, (alphabetsize+1)*sizeof(EnerGradX));
    int count = 0;
    for (int i = 0; i <= MAXATOMTYPES; i++) {
      if (i == MAXATOMTYPES || alphabet[i]) {
        p[i] = nr_energrads + count + 1;
        count++;
      }        
    }
    nr_energrads += alphabetsize + 1;    
  }  
  
}

void calculate_line(
  float gridspacing,
  int gridextension,
  int gridx,
  int gridy,
  int gridz,
  int nratoms,
  const Coor *xb,
  double neighbourdissq, 
  bool cdie, double ffelec, double plateaudissq,
  bool calc_pot,
  
  int z,
  Coor *dis,
  
  double *wer,  double *charges, int *atomtypes,
  const Parameters &rc, const Parameters &ac, 
  const Parameters &emin, const Parameters &rmin2, const iParameters &ipon,   
  int potshape, float swi_on, float swi_off,
  
  int alphabetsize,
  bool (&alphabet)[MAXATOMTYPES],
  
  Voxel *cinnergrid,
  Potential *cbiggrid,
  int &potcount,
  Neighbour *neighbours,
  int &nr_neighbours,
  EnerGradX *&energrads,
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
  if (z < 0 || z >= gridz) junction_z = 0;

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

  Neighbour *neigh = new Neighbour[NEIGHBOUR_MAX];  
  memset(neigh, 0, NEIGHBOUR_MAX*sizeof(Neighbour));
  
  for (int y = -gridextension; y < gridy+gridextension; y++) {
    inner_zy = inner_z;    
    if (y < 0 || y >= gridy) inner_zy = 0;
    junction_zy = junction_z;
    if (y < 0 || y >= gridy) junction_zy = 0;

    memcpy(diszyx, diszy, nratoms*sizeof(Coor));    
    memcpy(disjzyx, disjzy, nratoms*sizeof(Coor));    
    for (int x = -gridextension; x < gridx+gridextension; x++) {
      inner_zyx = inner_zy;    
      if (x < 0 || x >= gridx) inner_zyx = 0;
      junction_zyx = junction_zy;
      if (x < 0 || x >= gridx) junction_zyx = 0;
      if (inner_zyx) { 
        int index = x + gridx * y;
        Voxel &v = cinnergrid[index];
        if (calc_pot) {
          #pragma omp atomic
          potcount++;
          _calc_potential(v.potential, diszyx, xb, nratoms, 
           wer,charges,atomtypes,
           rc,ac,emin,rmin2,ipon,potshape,swi_on, swi_off,
           alphabetsize,alphabet,nr_energrads,energrads,cdie,ffelec,plateaudissq
          ); 
        }
        if (junction_zyx) {
          v.nr_neighbours = _calc_neighbours(
           v.neighbourlist, disjzyx, nratoms,atomtypes,
           neighbourdissq,plateaudissq,nr_neighbours, neighbours, neigh
          );
        }    
      }
      if (!(x % 2) && !(y % 2) && !(z % 2))
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
           alphabetsize,alphabet,nr_energrads,energrads,cdie,ffelec,plateaudissq
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
     gridspacing, \
     gridextension, \
     gridx, \
     gridy, \
     gridz, \
     nratoms, \
     xb, \
     neighbourdissq, \
     cdie, ffelec, plateaudissq, \
     calc_pot, \
     \
     z, \
     dis, \
     \
     wer,  charges, atomtypes, \
     rc, ac,  \
     emin, rmin2, ipon, \
     potshape, swi_on, swi_off, \
     alphabetsize, \
     alphabet, \
     \
     cinnergrid, \
     cbiggrid, \
     potcount, \
     neighbours, \
     nr_neighbours, \
     energradsX, \
     nr_energrads \
    ); 


inline void _calc_potential_elec(float &energy, Coorf &grad, 
#ifdef TORQUEGRID
float (&torques)[9],
#endif
const  Coor *dis, const Coor *xb, int nrdis, double *charges,
float swi_on, float swi_off, bool cdie, double ffelec, double plateaudissq

) {
  double energy0 = 0; Coor grad0 = {0,0,0};
  for (int n = 0; n < nrdis; n++) {
 
    float charge = charges[n];
    if (fabs(charge) > 0.001) {    
      const Coor &d = dis[n];
      double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
      if (dsq > distcutoffsq_elec) continue; 

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

      double rr2 = 1.0/dsq;
      double ratio = 1;
      if (dsq < plateaudissq) {
        ratio = sqrt(dsq/plateaudissq);
        dsq = plateaudissq;
        rr2 = 1.0/dsq;
      }
      Coor dd = {d[0]*rr2,d[1]*rr2,d[2]*rr2};          
      
      double charge0 = charge * ffelec * ffelec; //at realtime: to be multiplied with the charge of the other atom...
      
      double energy00; Coor grad00;
      elec(1,cdie,
        charge0,rr2,dd[0]*ratio,dd[1]*ratio,dd[2]*ratio,fswi,energy00,grad00);
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

inline void _calc_potential(
 Potential &p, const Coor *dis, const Coor *xb, int nrdis,
 double *wer, double *charges,int *atomtypes,
 const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
 const iParameters &ipon, const int &potshape, const float &swi_on, const float &swi_off, //ATTRACT params
 int alphabetsize,
 bool (&alphabet)[MAXATOMTYPES],
 int &nr_energrads, EnerGradX *&energrads,
 bool cdie, double ffelec, double plateaudissq //grid params
) {
  EnerGradX energrads0[MAXATOMTYPES+1];
  memset(energrads0, 0, sizeof(energrads0));
  
  for (int n = 0; n < nrdis; n++) {
    int atomtype = atomtypes[n];
    if (atomtype == 0) continue;  
    const Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2]; 
    if (dsq >= distcutoffsq_vdw) continue;

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
      //EnerGradX &e = *(energrads + p[i] - 1);
      EnerGradX &e = energrads0[i];
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
  //EnerGradX &e = *(energrads + p[MAXATOMTYPES] - 1);
  EnerGradX &e = energrads0[MAXATOMTYPES];
  _calc_potential_elec(e.energy, e.grad, 
#ifdef TORQUEGRID    
   e.torques,
#endif
   dis, xb, nrdis, charges, swi_on, swi_off, cdie, ffelec, plateaudissq   
  );  
  
  memset(&p, 0, sizeof(Potential));
  new_energrads(p, alphabetsize, alphabet, energrads0, nr_energrads, energrads);  
}

inline short _calc_neighbours(int &neighbourlist, const Coor *dis, int nrdis, int *atomtypes, 
double neighbourdissq, double plateaudissq,
int &nr_neighbours, Neighbour *neighbours, Neighbour *neigh
) {
  short neighboursize = 0;
  for (int n = 0; n < nrdis; n++) {
    if (atomtypes[n] == 0) continue;  
    const Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];       
    if (dsq >= neighbourdissq) continue;
    if (neighboursize == NEIGHBOUR_MAX) {fprintf(stderr, "Neighbour overflow\n"); exit(1);}
    Neighbour &nb = neigh[neighboursize];
    neighboursize++;    
    nb.type = 2;
    if (dsq <= plateaudissq) nb.type = 1;    
    nb.index = n;
  }

  if (neighboursize) {
    #pragma omp critical(neighbourlist)
    {  
      neighbourlist = nr_neighbours+1;
      if (nr_neighbours+neighboursize > MAXGRIDNEIGHBOUR) {fprintf(stderr, "Total neighbour grid size overflow\n"); exit(1);}
      memcpy(neighbours+nr_neighbours, neigh, neighboursize*sizeof(Neighbour));
      nr_neighbours += neighboursize;
    }
  }
  
  return neighboursize;
}
#ifdef TORQUEGRID
void Grid::calculate_torque(
#else
void Grid::calculate_std(
#endif  
  int cartstatehandle, double plateaudis, double neighbourdis, float gridspacing, int gridextension, bool (&alphabet)[MAXATOMTYPES], bool cdie, float epsilon, bool calc_pot) {  
  init(gridspacing, gridextension, plateaudis, neighbourdis, alphabet);  

  double ffelec = sqrt(felec/epsilon);
  
  nr_neighbours = 0; 
  shm_neighbours = -1; 
  energradsX = NULL;
  nr_energrads = 0;
  shm_energrads = -1;
  
  int n;
  
  
  //Load Cartesian state and parameters
  CartState &cartstate = cartstate_get(cartstatehandle); 
  int nratoms = cartstate.natom[0];
  int start = 0;
  Coor *x = (Coor *) &(cartstate.x[0]);
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
      
  //Calculate grid dimensions
  double minx=99999, miny=99999, minz=99999;
  double maxx=-99999, maxy=-99999, maxz=-99999;
  for (n = 0; n < nratoms; n++) {
    Coor &a = x[n];
    if (a[0] > maxx) maxx = a[0]; if (a[0] < minx) minx = a[0];
    if (a[1] > maxy) maxy = a[1]; if (a[1] < miny) miny = a[1];
    if (a[2] > maxz) maxz = a[2]; if (a[2] < minz) minz = a[2];    
  }
  
    
  //Store some parameters
  this->torquegrid = IS_TORQUEGRID;
  this->gridx = ceil((maxx-minx+2*plateaudis)/gridspacing);
  this->gridy = ceil((maxy-miny+2*plateaudis)/gridspacing);
  this->gridz = ceil((maxz-minz+2*plateaudis)/gridspacing);    
  this->pivot[0] = cartstate.pivot[0][0];
  this->pivot[1] = cartstate.pivot[1][0];
  this->pivot[2] = cartstate.pivot[2][0];  
  this->natoms = nratoms;
  this->ori[0] = minx-plateaudis;
  this->ori[1] = miny-plateaudis;
  this->ori[2] = minz-plateaudis;
  
  //Initialize grids
  gridx2 = int((gridx + 2 * gridextension)/2)+1;
  gridy2 = int((gridy + 2 * gridextension)/2)+1;
  gridz2 = int((gridz + 2 * gridextension)/2)+1;

  Voxel *cinnergrid;
  Potential *cbiggrid;
  
  this->innergrid = new Voxel[gridx*gridy*gridz];
  memset(innergrid, 0, gridx*gridy*gridz*sizeof(Voxel));

  if (calc_pot) {
    this->biggrid = new Potential[gridx2*gridy2*gridz2];
    memset(biggrid, 0, gridx2*gridy2*gridz2*sizeof(Potential));
  }
  this->neighbours = new Neighbour[MAXGRIDNEIGHBOUR];
  if (neighbours == NULL) {
    fprintf(stderr, "Could not allocate memory for %d Neighbours\n", MAXGRIDNEIGHBOUR);
    exit(1);
  }
  memset(this->neighbours, 0,MAXGRIDNEIGHBOUR*sizeof(Neighbour));
  
  //Set up distances 
  Coor *dis = new Coor[nratoms]; 
  for (n = 0; n < nratoms; n++) {
    dis[n][0] = -(x[n][0]-ori[0]);
    dis[n][1] = -(x[n][1]-ori[1]);
    dis[n][2] = -(x[n][2]-ori[2]);
  }          

  //Main loop
  int potcount = 0;
  #pragma omp parallel for private(cinnergrid,cbiggrid) schedule(dynamic,1) 
  for (int z = -gridextension; z < gridz+gridextension; z++) {    
    fprintf(stderr, "%d/%d %d %d\n", z, gridz+gridextension, potcount, nr_neighbours);  
        
    cinnergrid=NULL;
    if (z >= 0) cinnergrid = innergrid + gridx*gridy*z;      
    int zz = (z + gridextension)/2;
    cbiggrid = biggrid + gridx2*gridy2*zz;
    
    MACRO_calculate_line    
  }     
  fprintf(stderr, "%d/%d %d %d\n", gridz+gridextension, gridz+gridextension, potcount, nr_neighbours);    

}
