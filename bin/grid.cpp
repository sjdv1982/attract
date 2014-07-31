#include "grid.h"
#include "nonbon.h"
#include "state.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <unistd.h>

void error(const char *filename) {
  fprintf(stderr, "Reading error in grid file %s\n", filename);
  exit(1);
}

void get_shm_name(int shm_id, char *shm_name) {
  sprintf(shm_name, "/attract-grid%d", shm_id);
}

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
    _ratio[n] = sqrt(dissq/plateaudissq) / (1/dissq) * (1/plateaudissq);
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

void Grid::read(const char *filename) {
  int n;  

  int read;
  
  FILE *f = fopen(filename, "rb");  
  if (f == NULL) error(filename);

  read = fread(&architecture, sizeof(short),1,f);  
  if (architecture != ARCHITECTURE) {
    fprintf(stderr, "Reading error in grid file %s, grid was computed on %d bit, but we are on %d bit\n", filename, architecture, ARCHITECTURE);
    exit(1);
  }  
  if (!read) error(filename);  
  read = fread(&gridspacing, sizeof(double),1,f);  
  if (!read) error(filename);
  read = fread(&gridextension, sizeof(int),1,f);
  if (!read) error(filename);
  read = fread(&plateaudis, sizeof(double),1,f);  
  if (!read) error(filename);
  read = fread(&neighbourdis, sizeof(double),1,f);    
  if (!read) error(filename);
  bool alphabet0[MAXATOMTYPES];
  read = fread(&alphabet0, sizeof(alphabet0),1,f);    
  if (!read) error(filename);
  init(gridspacing, gridextension, plateaudis,neighbourdis,alphabet0);

  float arr1[3];
  read = fread(arr1, 3*sizeof(float),1,f);
  if (!read) error(filename);
  ori[0]=arr1[0];ori[1]=arr1[1];ori[2]=arr1[2];
  int arr2[6];
  read = fread(arr2, 6*sizeof(int),1,f);
  if (!read) error(filename);
  gridx=arr2[0];gridy=arr2[1];gridz=arr2[2];  
  gridx2=arr2[3];gridy2=arr2[4];gridz2=arr2[5];  
  read = fread(&natoms,sizeof(int),1,f);
  if (!read) error(filename);
  read=fread(&pivot,sizeof(Coor),1,f);
  if (!read) error(filename);
  
  read = fread(&nr_energrads, sizeof(nr_energrads),1,f);  
  if (!read) error(filename);
  read = fread(&shm_energrads, sizeof(shm_energrads),1,f);  
  if (!read) error(filename);

  if (nr_energrads) {
    if (shm_energrads == -1) {
      energrads = new EnerGrad[nr_energrads];
      read = fread(energrads, nr_energrads*sizeof(EnerGrad),1,f);  
      if (!read) error(filename);
    }
    else {
      char shm_name[100];
      get_shm_name(shm_energrads, shm_name);  
      int fshm1 = shm_open(shm_name, O_RDONLY, S_IREAD);
      if (fshm1 == -1) {
	fprintf(stderr, "Reading error in grid file %s: shared memory segment %d for potential list does not exist\n", filename, shm_energrads);
	exit(1); 
      }   
      ftruncate(fshm1, nr_energrads*sizeof(EnerGrad));
      energrads = (EnerGrad *) mmap(0,nr_energrads*sizeof(EnerGrad),
       PROT_READ, MAP_SHARED | MAP_NORESERVE, fshm1, 0);
      if (energrads == NULL) {
	fprintf(stderr, "Reading error in grid file %s: Could not load shared memory segment %d\n", filename, shm_energrads);
	exit(1);       
      }  
    }  
  }
  
  read = fread(&nr_neighbours, sizeof(nr_neighbours),1,f);  
  if (!read) error(filename);
  read = fread(&shm_neighbours, sizeof(shm_neighbours),1,f);  
  if (!read) error(filename);
  if (shm_neighbours == -1) {
    neighbours = new Neighbour[nr_neighbours];
    read = fread(neighbours, nr_neighbours*sizeof(Neighbour),1,f);
    if (!read) error(filename);
  }
  else {
    char shm_name[100];
    get_shm_name(shm_neighbours, shm_name);  
    int fshm2 = shm_open(shm_name, O_RDONLY, S_IREAD);
    if (fshm2 == -1) {
      fprintf(stderr, "Reading error in grid file %s: shared memory segment %d for neighbour list does not exist\n", filename, shm_neighbours);
      exit(1); 
    }   
    ftruncate(fshm2, nr_neighbours*sizeof(Neighbour));
    neighbours = (Neighbour *) mmap(0,nr_neighbours*sizeof(Neighbour),
     PROT_READ, MAP_SHARED | MAP_NORESERVE, fshm2, 0);
    if (neighbours == NULL) {
      fprintf(stderr, "Reading error in grid file %s: Could not load shared memory segment %d\n", filename, shm_neighbours);
      exit(1); 
    }  
  }  
  long innergridsize, biggridsize;
  read = fread(&innergridsize, sizeof(innergridsize),1,f);
  if (!read) error(filename);
  innergrid = new Voxel[innergridsize];
  read = fread(innergrid, innergridsize*sizeof(Voxel),1,f);
  if (!read) error(filename);
  read = fread(&biggridsize, sizeof(biggridsize),1,f);
  if (!read) error(filename);
  if (biggridsize) {
    biggrid = new Potential[biggridsize];
    read = fread(biggrid, biggridsize*sizeof(Potential),1,f);
  }
  if (!read) error(filename);
  fclose(f);
  for (n = 0; n < innergridsize; n++) {
    if (innergrid[n].potential[MAXATOMTYPES]) {
      for (int i = 0; i <= MAXATOMTYPES; i++) {
        if (i < MAXATOMTYPES && alphabet[i] == 0) continue;
        int dif = innergrid[n].potential[i] - 1;
	if (dif < 0 || dif  >= nr_energrads) {
	  fprintf(stderr, "Reading error in %s, innergrid voxel %d atom type %d: %d >= %d\n",             
	  filename, n+1, i+1, dif, nr_energrads);
	  exit(1);
	}	
      }
      int nl = innergrid[n].neighbourlist;
      int nr = innergrid[n].nr_neighbours;
      bool empty1 = (nl == 0);
      bool empty2 = (nr == 0);
      
      if (empty1 != empty2 || (nl + nr - 1 > nr_neighbours)) {
	fprintf(stderr, "Reading error in %s, innergrid voxel %d neighbourlist: %d + %d >= %d\n",             
	filename, n+1, nl-1, nr, nr_neighbours);
	exit(1);
      }	

    }  
  }
  for (n = 0; n < biggridsize; n++) {
    if (biggrid[n][MAXATOMTYPES]) {
      for (int i = 0; i <= MAXATOMTYPES; i++) {
        if (i < MAXATOMTYPES && alphabet[i] == 0) continue;
        int dif = biggrid[n][i]-1;
	if (dif < 0 || dif >= nr_energrads) {
	  fprintf(stderr, "Reading error in %s, biggrid voxel %d atom type %d: %d >= %d\n",             
	  filename, n+1, i+1, dif, nr_energrads);
	  exit(1);
	}
      }
    }
  }

  init(gridspacing, gridextension, plateaudis,neighbourdis,alphabet0);
}

void Grid::write(const char *filename) {
  long n;
  long innergridsize = gridx*gridy*gridz;
  long biggridsize = gridx2*gridy2*gridz2;    
  if (!nr_energrads) biggridsize = 0;
 
  FILE *f = fopen(filename, "wb");
  if (f == NULL) {
    fprintf(stderr, "Grid::write error for %s: Cannot open file for writing\n", filename);
    exit(1);
  }

  
  EnerGrad *shmptr1 = NULL; Neighbour *shmptr2 = NULL;
  if (nr_energrads && shm_energrads != -1) {
    char shm_name[100];
    get_shm_name(shm_energrads, shm_name);
    int fshm1 = shm_open(shm_name, (O_CREAT | O_RDWR), (S_IREAD | S_IWRITE));
    if (fshm1 == -1) {
      fprintf(stderr, "Grid::write error for %s: Cannot open shared memory for writing\n", filename);
      exit(1);
    }     
    ftruncate(fshm1, nr_energrads*sizeof(EnerGrad));
    shmptr1 = (EnerGrad *) mmap(0,nr_energrads*sizeof(EnerGrad),(PROT_READ | PROT_WRITE),
     MAP_SHARED , fshm1, 0);
    if (shmptr1 == MAP_FAILED) {
      fprintf(stderr, "Grid::write error for %s: Cannot map shared memory for writing\n", filename);
      exit(1);
    } 
    memset(shmptr1,0,nr_energrads*sizeof(EnerGrad));
    close(fshm1);    
  }
  
  if (shm_neighbours != -1) {
    char shm_name[100];
    get_shm_name(shm_neighbours, shm_name);
    int fshm2 = shm_open(shm_name, (O_CREAT | O_RDWR), (S_IREAD | S_IWRITE));
    if (fshm2 == -1) {
      fprintf(stderr, "Grid::write error for %s: Cannot open shared memory for writing\n", filename);
      exit(1);
    }     
    ftruncate(fshm2, nr_neighbours*sizeof(Neighbour));
    shmptr2 = (Neighbour *) mmap(0,nr_neighbours*sizeof(Neighbour),(PROT_READ | PROT_WRITE),
     MAP_SHARED, fshm2, 0);
    if (shmptr2 == MAP_FAILED) {
      fprintf(stderr, "Grid::write error for %s: Cannot map shared memory for writing\n", filename);
      exit(1);
    } 
    memset(shmptr2,0,nr_neighbours*sizeof(Neighbour));
    close(fshm2);
  }  
    
  fwrite(&architecture, sizeof(unsigned short),1,f);      
  fwrite(&gridspacing, sizeof(double),1,f);  
  fwrite(&gridextension, sizeof(int),1,f);
  fwrite(&plateaudis, sizeof(double),1,f);  
  fwrite(&neighbourdis, sizeof(double),1,f);    
  fwrite(&alphabet, sizeof(alphabet),1,f);
  float arr1[] = {ori[0],ori[1],ori[2]};
  fwrite(arr1, 3*sizeof(float),1,f);
  int arr2[] = {gridx,gridy,gridz,gridx2,gridy2,gridz2};
  fwrite(arr2, 6*sizeof(int),1,f);
  fwrite(&natoms, sizeof(int),1,f); 
  fwrite(&pivot,sizeof(Coor),1,f);

  if (nr_energrads) {
    energrads = (EnerGrad *) realloc(energrads, nr_energrads*sizeof(EnerGrad));
    EnerGrad *energrads_reordered = new EnerGrad[nr_energrads];
    memset(energrads_reordered, 0, nr_energrads*sizeof(EnerGrad));
    if (energrads_reordered) { //only re-order the energrads if we got the memory for it    
      int nr_energrads2 = 0;  
      for (n = 0; n < innergridsize; n++) {
        if (innergrid[n].potential[MAXATOMTYPES]) {
          for (int i = 0; i <= MAXATOMTYPES; i++) {
            if (i < MAXATOMTYPES && alphabet[i] == 0) continue;
            unsigned int &oldpos = innergrid[n].potential[i];
            unsigned int newpos = nr_energrads2+1;
            memcpy(&energrads_reordered[newpos-1], &energrads[oldpos-1], sizeof(EnerGrad));
            oldpos = newpos;
            nr_energrads2++;           
          }
        }  
      }
      for (n = 0; n < biggridsize; n++) {
        if (biggrid[n][MAXATOMTYPES]) {
          for (int i = 0; i <= MAXATOMTYPES; i++) {
            if (i < MAXATOMTYPES && alphabet[i] == 0) continue;
            unsigned int &oldpos = biggrid[n][i];
            unsigned int newpos = nr_energrads2+1;
            memcpy(&energrads_reordered[newpos-1], &energrads[oldpos-1], sizeof(EnerGrad));
            oldpos = newpos;
            nr_energrads2++;           
          }
        }  
      }  
      if (nr_energrads != nr_energrads2) {
        fprintf(stderr, "ERR nr_energrads %d %d\n", nr_energrads, nr_energrads2);
      }
      free(energrads);
      energrads = energrads_reordered;      
    }  
  }
  
  fwrite(&nr_energrads, sizeof(nr_energrads),1,f);
  fwrite(&shm_energrads, sizeof(shm_energrads),1,f);
   
  if (nr_energrads) { 
    if (shm_energrads == -1) 
      fwrite(energrads,nr_energrads*sizeof(EnerGrad),1,f);
    else 
      memcpy(shmptr1, energrads,nr_energrads*sizeof(EnerGrad));
  }  
  
  if (nr_neighbours) {
    Neighbour *neighbours_reordered = new Neighbour[nr_neighbours];
    memset(neighbours_reordered, 0, nr_neighbours*sizeof(Neighbour));
    if (neighbours_reordered) { //only re-order the neighbours if we got the memory for it    
      int nr_neighbours2 = 0;  
      for (n = 0; n < innergridsize; n++) {
        Voxel &v = innergrid[n];
        if (v.nr_neighbours) {
          memcpy(&neighbours_reordered[nr_neighbours2], &neighbours[v.neighbourlist-1], v.nr_neighbours*sizeof(Neighbour));
          v.neighbourlist = nr_neighbours2+1;
          nr_neighbours2 += v.nr_neighbours;
        }        
      }
      if (nr_neighbours2 != nr_neighbours) {
        fprintf(stderr, "ERR nr_neighbours %d %d\n", nr_neighbours, nr_neighbours2);
      }      
      delete[] neighbours;
      neighbours = neighbours_reordered;                
    }
  }
  
  fwrite(&nr_neighbours, sizeof(nr_neighbours),1,f);
  fwrite(&shm_neighbours, sizeof(shm_neighbours),1,f);
  if (shm_neighbours == -1) 
    fwrite(neighbours, nr_neighbours*sizeof(Neighbour), 1, f);
  else 
    memcpy(shmptr2, neighbours, nr_neighbours*sizeof(Neighbour));
    
  fwrite(&innergridsize, sizeof(innergridsize),1,f);
  fwrite(innergrid, innergridsize*sizeof(Voxel),1,f);
  fwrite(&biggridsize, sizeof(biggridsize),1,f);
  if (biggridsize) {
    fwrite(biggrid, biggridsize*sizeof(Potential),1,f);
  }
  fclose(f);
}

#ifdef TORQUEGRID
inline void add_potential(EnerGrad *energrads, Potential &p, int iab, double charge, int atomtype, double w, int fixre, double &evdw, double &eelec, Coor &grad, double *deltar, double (&rtorques)[3][3]) {
#else
inline void add_potential(EnerGrad *energrads, Potential &p, int iab, double charge, int atomtype, double w, int fixre, double &evdw, double &eelec, Coor &grad, double *deltar) {
#endif
  if (p[atomtype] == 0) {
    evdw  += w * 999999;
  }
  else {
    EnerGrad &e = *(energrads+p[atomtype]-1);
    evdw  += w * e.energy;
    if (iab) {
      Coor dgrad = {
        w * e.grad[0],
	w * e.grad[1],
	w * e.grad[2]
      };
      grad[0] += dgrad[0];
      grad[1] += dgrad[1];
      grad[2] += dgrad[2];
      if (!fixre) {
	deltar[3] -= dgrad[0];
	deltar[4] -= dgrad[1];
	deltar[5] -= dgrad[2];
#ifdef TORQUEGRID
	rtorques[0][0] += w * e.torques[0];
	rtorques[0][1] += w * e.torques[1];
	rtorques[0][2] += w * e.torques[2];
	rtorques[1][0] += w * e.torques[3];
	rtorques[1][1] += w * e.torques[4];
	rtorques[1][2] += w * e.torques[5];
	rtorques[2][0] += w * e.torques[6];
	rtorques[2][1] += w * e.torques[7];
	rtorques[2][2] += w * e.torques[8];
#endif      
      }
    }
    if (fabs(charge) > 0.001) {
      double wcharge = w * charge;
      EnerGrad &e = *(energrads+p[MAXATOMTYPES]-1);
      eelec  += wcharge * e.energy;
      if (iab) {
	Coor dgrad = {
          wcharge * e.grad[0],
	  wcharge * e.grad[1],
	  wcharge * e.grad[2]
	};
	grad[0] += dgrad[0];
	grad[1] += dgrad[1];
	grad[2] += dgrad[2];
	if (!fixre) {
	  deltar[3] -= dgrad[0];
	  deltar[4] -= dgrad[1];
	  deltar[5] -= dgrad[2];
#ifdef TORQUEGRID	
	  rtorques[0][0] += wcharge * e.torques[0];
	  rtorques[0][1] += wcharge * e.torques[1];
	  rtorques[0][2] += wcharge * e.torques[2];
	  rtorques[1][0] += wcharge * e.torques[3];
	  rtorques[1][1] += wcharge * e.torques[4];
	  rtorques[1][2] += wcharge * e.torques[5];
	  rtorques[2][0] += wcharge * e.torques[6];
	  rtorques[2][1] += wcharge * e.torques[7];
	  rtorques[2][2] += wcharge * e.torques[8];      	
#endif	
        }
      }	
    }
  }
}

#ifdef TORQUEGRID
inline void add_potential2(EnerGrad *energrads, Potential p, int iab, double charge, int atomtype, double w, int fixre, double &evdw, double &eelec, Coor &grad, double *deltar, double (&rtorques)[3][3]) {
#else
inline void add_potential2(EnerGrad *energrads, Potential p, int iab, double charge, int atomtype, double w, int fixre, double &evdw, double &eelec, Coor &grad, double *deltar) {
#endif
  if (p[atomtype]) {
    EnerGrad &e = *(energrads + p[atomtype] - 1);
    evdw  += w * e.energy;
    if (iab) {
      Coor dgrad = {
        w * e.grad[0],
	w * e.grad[1],
	w * e.grad[2]
      };
      grad[0] += dgrad[0];
      grad[1] += dgrad[1];
      grad[2] += dgrad[2];
      if (!fixre) {
	deltar[3] -= dgrad[0];
	deltar[4] -= dgrad[1];
	deltar[5] -= dgrad[2];
#ifdef TORQUEGRID      
	rtorques[0][0] += w * e.torques[0];
	rtorques[0][1] += w * e.torques[1];
	rtorques[0][2] += w * e.torques[2];
	rtorques[1][0] += w * e.torques[3];
	rtorques[1][1] += w * e.torques[4];
	rtorques[1][2] += w * e.torques[5];
	rtorques[2][0] += w * e.torques[6];
	rtorques[2][1] += w * e.torques[7];
	rtorques[2][2] += w * e.torques[8];      
#endif      
     }
    }
    if (fabs(charge) > 0.001) {
      double wcharge = w * charge;
      EnerGrad &e = *(energrads + p[MAXATOMTYPES] - 1);
      eelec  += wcharge * e.energy;
      if (iab) {
	Coor dgrad = {
          wcharge * e.grad[0],
	  wcharge * e.grad[1],
	  wcharge * e.grad[2]
	};
	grad[0] += dgrad[0];
	grad[1] += dgrad[1];
	grad[2] += dgrad[2];
	if (!fixre) {	
	  deltar[3] -= dgrad[0];
	  deltar[4] -= dgrad[1];
	  deltar[5] -= dgrad[2];
#ifdef TORQUEGRID	
	  rtorques[0][0] += wcharge * e.torques[0];
	  rtorques[0][1] += wcharge * e.torques[1];
	  rtorques[0][2] += wcharge * e.torques[2];
	  rtorques[1][0] += wcharge * e.torques[3];
	  rtorques[1][1] += wcharge * e.torques[4];
	  rtorques[1][2] += wcharge * e.torques[5];
	  rtorques[2][0] += wcharge * e.torques[6];
	  rtorques[2][1] += wcharge * e.torques[7];
	  rtorques[2][2] += wcharge * e.torques[8];	
#endif	
        }
      }
    }    
  }
}

#ifdef TORQUEGRID
inline void trilin(
  const Coor &d, int atomtype, double charge0, double charge, double wel, //ligand
  const Grid &g, bool rigid, const Coor *xr, const double *wer, const double *chair, const int *iacir,  //receptor

  const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, int potshape, int cdie, float swi_on, float swi_off, 
  //ATTRACT params

  int iab, int fixre, double &evdw, double &eelec, Coor &grad, //output

  Coor *fr, const double (&pm2)[3][3][3], double *deltar, double (&rtorques)[3][3] //receptor
) 
#define TORQUEARGS ,rtorques  
#else
inline void trilin(
  const Coor &d, int atomtype, double charge0, double charge, double wel, //ligand
  const Grid &g, bool rigid, const Coor *xr, const double *wer, const double *chair, const int *iacir,  //receptor

  const Parameters &rc, const Parameters &ac, const Parameters &emin, const Parameters &rmin2,
  const iParameters &ipon, int potshape, int cdie, float swi_on, float swi_off,  //ATTRACT params

  int iab, int fixre, double &evdw, double &eelec, Coor &grad, //output
  
  Coor *fr, double *deltar //receptor
) 
#define TORQUEARGS 
#endif
{
  bool has_pot = (g.nr_energrads > 0);
  
  double ax = (d[0]-g.ori[0])/g.gridspacing;
  double ay = (d[1]-g.ori[1])/g.gridspacing;
  double az = (d[2]-g.ori[2])/g.gridspacing;
  //std::cerr << "Ligand position " << d[0] << " " << d[1] << " " << d[2] << "\n";
  //std::cerr << "Grid " << ax << " " << ay << " " << az << "\n";
  bool chargenonzero = (fabs(charge0) > 0.001);

  double aax = (ax+g.gridextension)/2;  
  if (aax < 0 || aax >= g.gridx2-1) return;
  double aay = (ay+g.gridextension)/2;
  if (aay < 0 || aay >= g.gridy2-1) return;
  double aaz = (az+g.gridextension)/2;  
  if (aaz < 0 || aaz >= g.gridz2-1) return;
  
  bool inside = ( 
    (ax >= 0 && ax < g.gridx-1) &&
    (ay >= 0 && ay < g.gridy-1) &&
    (az >= 0 && az < g.gridz-1) 
  );
  int gx,gy;
  if (!inside) {
    ax = aax; ay = aay; az = aaz;
    gx = g.gridx2; gy = g.gridy2; 
  }
  else {
    gx = g.gridx; gy = g.gridy; 
  }
 // std::cerr << "Grid points " << ax << " " << ay << " " << az << "\n";
  if (g.nr_energrads) {
    int px0 = floor(ax);  
    int px1 = ceil(ax);
    double wx1 = ax - px0;

    int py0 = floor(ay);  
    int py1 = ceil(ay);
    double wy1 = ay - py0;

    int pz0 = floor(az);  
    int pz1 = ceil(az);
    double wz1 = az - pz0;    
    double wx0 = 1-wx1, wy0 = 1-wy1, wz0 = 1-wz1;

    wx0 *= wel; wy0 *= wel; wz0 *= wel;
    wx1 *= wel; wy1 *= wel; wz1 *= wel;

    Potential *p;
 //   std::cerr << "Add potential px0 etc " << px0 << " " << px1 << " " <<py0 << " " << py1 << " " << pz0 << " " << pz1 << "\n";
    if (inside) {
      p = &g.innergrid[px0+gx*py0+gx*gy*pz0].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx0*wy0*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
   //   std::cerr << evdw << " " << px0+gx*py0+gx*gy*pz0 << "\t";
      p = &g.innergrid[px0+gx*py0+gx*gy*pz1].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx0*wy0*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
  //    std::cerr << evdw << " " << px0+gx*py0+gx*gy*pz1<< "\t";
      p = &g.innergrid[px0+gx*py1+gx*gy*pz0].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx0*wy1*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
  //    std::cerr << evdw << " " << px0+gx*py1+gx*gy*pz0<< "\t";
      p = &g.innergrid[px0+gx*py1+gx*gy*pz1].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx0*wy1*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
 //     std::cerr << evdw << " " << px0+gx*py1+gx*gy*pz1<< "\t";
      p = &g.innergrid[px1+gx*py0+gx*gy*pz0].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx1*wy0*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
  //    std::cerr << evdw << "\t";
      p = &g.innergrid[px1+gx*py0+gx*gy*pz1].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx1*wy0*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
//      std::cerr << evdw << "\t";
      p = &g.innergrid[px1+gx*py1+gx*gy*pz0].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx1*wy1*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
  //    std::cerr << evdw << "\t";
      p = &g.innergrid[px1+gx*py1+gx*gy*pz1].potential;
      add_potential(g.energrads, *p,iab,charge,atomtype,wx1*wy1*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
    }
    else {
      p = &g.biggrid[px0+gx*py0+gx*gy*pz0];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx0*wy0*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);    
      p = &g.biggrid[px0+gx*py0+gx*gy*pz1];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx0*wy0*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px0+gx*py1+gx*gy*pz0];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx0*wy1*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px0+gx*py1+gx*gy*pz1];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx0*wy1*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px1+gx*py0+gx*gy*pz0];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx1*wy0*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px1+gx*py0+gx*gy*pz1];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx1*wy0*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px1+gx*py1+gx*gy*pz0];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx1*wy1*wz0, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
      p = &g.biggrid[px1+gx*py1+gx*gy*pz1];
      add_potential2(g.energrads, *p,iab,charge,atomtype,wx1*wy1*wz1, fixre, evdw, eelec, grad ,deltar TORQUEARGS);
    }
  }
 // std::cerr << "After adding potential " << evdw << "\n";
  if (inside) {
    int pxj = floor(ax+0.5);
    int pyj = floor(ay+0.5);
    int pzj = floor(az+0.5);
    if (pxj < g.gridx && pyj < g.gridy && pzj < g.gridz) {
      int index = pxj+gx*pyj+gx*gy*pzj;
      Voxel &v = g.innergrid[index];
      if (v.nr_neighbours > 0) {
        for (int i = 0; i < v.nr_neighbours; i++) {
	  Neighbour &nb = g.neighbours[v.neighbourlist+i-1];
	  if (rigid && nb.type == 2) continue;
	  const Coor &dd = xr[nb.index];

	  Coor dis = {d[0] - dd[0], 
	              d[1] - dd[1], 
		      d[2] - dd[2], 
		      };
	  
	  double dsq = dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
	  if (dsq < 2.0){
	//   std::cerr << nb.index << "\t" << d[0] << "\t" << d[1] << "\t" << d[2] << "\t" << dd[0] << "\t" << dd[1] << "\t" << dd[2] << "\n"; 
	  }
          if (dsq > g.plateaudissq) continue;
	  int atomtype2 = iacir[nb.index]-1;
	  double ww = wel * wer[nb.index];

  	  Coor &gradr = fr[nb.index];
          double fswi = 1, fswi2 = 1;
          if (chargenonzero) {
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
	      if (g.plateaudissq > swi_on*swi_on) {
		if (g.plateaudissq > swi_off*swi_off) {
		  fswi2 = 0;
		}
		else {
		  double distance = sqrt(g.plateaudissq) ;
		  fswi2 = 1-(distance - swi_on)/(swi_off-swi_on);
		}
	      }    
                 
	    }
          }
          double rci = rc[atomtype][atomtype2];
          double aci = ac[atomtype][atomtype2];
          double emini = emin[atomtype][atomtype2];
          double rmin2i = rmin2[atomtype][atomtype2];
          int ivor = ipon[atomtype][atomtype2];
        
        
          double rr2 = 1.0/dsq;	    
          dis[0] *= rr2; dis[1] *= rr2; dis[2] *= rr2;

          double evdw0; Coor grad0;
          nonbon(iab,ww,rci,aci,emini,rmin2i,ivor, dsq, rr2, 
            dis[0], dis[1], dis[2], potshape, fswi, evdw0, grad0);	
          evdw += evdw0;
          grad[0] += grad0[0];
          grad[1] += grad0[1];
          grad[2] += grad0[2];
          if (!fixre) {
            gradr[0] -= grad0[0];
            gradr[1] -= grad0[1];
            gradr[2] -= grad0[2];
          }
//	    std::cerr << "New vdW energy " << evdw << " = oldenergy + " << evdw0 << "\n";
          bool calc_elec = 0;
          double c;
          if (chargenonzero) {
            c = charge0 * chair[nb.index] * ww;	    
            if (fabs(c) > 0.001) calc_elec = 1;
          }
          
          Coor rdis;
          if (has_pot || calc_elec) {
            double r = g.get_ratio(dsq);
            rdis[0] = r*dis[0];
            rdis[1] = r*dis[1];
            rdis[2] = r*dis[2];              
          }
          if (has_pot) {
            nonbon(iab,ww,rci,aci,emini,rmin2i,ivor, g.plateaudissq, 
              g.plateaudissqinv,
              rdis[0], rdis[1], rdis[2], potshape, fswi2, evdw0, grad0);
            evdw -= evdw0;
            grad[0] -= grad0[0];
            grad[1] -= grad0[1];
            grad[2] -= grad0[2];
            if (!fixre) {
              gradr[0] += grad0[0];
              gradr[1] += grad0[1];
              gradr[2] += grad0[2];
            }
          }  
//	    std::cerr << "New vdW energy " << evdw << " = oldenergy - " << evdw0 << "\n";
          if (calc_elec) {
            double eelec00; Coor grad00;
            elec(iab,cdie,c,rr2,dis[0],dis[1],dis[2],fswi,eelec00,grad00);
            eelec += eelec00;
            grad[0] += grad00[0];
            grad[1] += grad00[1];
            grad[2] += grad00[2];
            if (!fixre) {	    
              gradr[0] -= grad00[0];
              gradr[1] -= grad00[1];
              gradr[2] -= grad00[2];
            }
            elec(iab,cdie,c,g.plateaudissqinv,
              rdis[0],rdis[1],rdis[2],fswi2,eelec00,grad00);
            eelec -= eelec00;
            grad[0] -= grad00[0];
            grad[1] -= grad00[1];
            grad[2] -= grad00[2];
            if (!fixre) {	    		
              gradr[0] += grad00[0];
              gradr[1] += grad00[1];
              gradr[2] += grad00[2];
            }            
          }
	  
	}
      }
    } 
    
  }//end inside
}

extern "C" void nonbon_grid_(
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
  
  Coor *fr, const double (&pm2)[3][3][3], double *deltar)
{
  
  double ffelec = sqrt(felec/epsilon);
  
  evdw = 0;
  eelec = 0;
#ifdef TORQUEGRID  
  double rtorques[3][3];
  memset(rtorques,0,3*3*sizeof(double));
#endif  
  int counter2 = 0;
  for (int n = 0; n < natoml; n++) {
    int atomtype = iacil[n]-1;
    if (atomtype == -1) continue;
    counter2++;
    double charge0 = chail[n];
    double charge = charge0/ffelec;
    Coor &grad = fl[n];
    const Coor &d = xl[n];
   // std::cerr << "Processing atoms " << n << " " << d[0] << " " << d[1] << " " << d[2] << " Receptor: ";
#ifdef TORQUEGRID    
    trilin(xl[n],atomtype,charge0,charge,wel[n],*g,rigid,xr,wer,chair,iacir,
      rc,ac,emin,rmin2,ipon,potshape,cdie,swi_on,swi_off, iab,fixre,evdw,eelec,grad,fr,pm2,deltar,rtorques);
#else      
    trilin(xl[n],atomtype,charge0,charge,wel[n],*g,rigid,xr,wer,chair,iacir,
      rc,ac,emin,rmin2,ipon,potshape,cdie,swi_on,swi_off,iab,fixre,evdw,eelec,grad,fr,deltar);
#endif
 //   std::cerr << "Atom energy" << n << " " << evdw << " " << eelec << "\n";
  }
 // std::cerr << "Energy from grid " << evdw << " " << eelec << "\n";
#ifdef TORQUEGRID      
  for (int j = 0; j < 3; j++) { //euler angle
    for (int k = 0; k < 3; k++) {//gradient
      for (int l = 0; l < 3; l++) { //coordinate
        deltar[j] += rtorques[l][k] * pm2[l][j][k];
      }
    }
  } 
#endif
}  
  
