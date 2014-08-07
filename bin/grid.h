#ifndef GRID_H 
#define GRID_H /* to make sure that we include only once... */

#include  "max.h"
#include  <cstdio>

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ARCHITECTURE 64
#else
#define ARCHITECTURE 32
#endif
#else
#define ARCHITECTURE 64 //if not GCC, by default assume 64 bit
#endif

typedef double Coor[3]; /*Coor is double[3]*/
typedef float Coorf[3]; /*Coorf is float[3]*/

struct EnerGradStd {
  float energy;
  Coorf grad;
};

struct EnerGradTorque {
  float energy;
  Coorf grad;
  float torques[9];
};

typedef unsigned int Potential[MAXATOMTYPES+1]; 
/*A Potential consists of 100 indices-to-EnerGrad:
 99 for the atom types, the last one is for elec*/

struct Neighbour {
  char type; /*1 = under plateaudis, 2 = between plateaudis and neighbourdis*/
  unsigned int index; /*receptor atom index*/
};

struct Voxel {
  Potential potential;
  int neighbourlist;
  short nr_neighbours;
};

struct Grid {
  bool torquegrid; //are we a torque grid?
  unsigned short architecture; //32 or 64 bits
  double gridspacing; /* 0.9 A */
  int gridextension; /* 32; always make this number even!*/
  double plateaudis; /* 10 A */
  double plateaudissq;
  double plateaudissqinv;
  double neighbourdis; 
  double neighbourdissq;
  int natoms;
  bool alphabet[MAXATOMTYPES];
  int alphabetsize; //the number of non-zero elements in alphabet

  Coor pivot;
  Coor ori; 
  /*origin of the inner grid (position of voxel (0,0,0) )*/
    
  int gridx, gridy, gridz; 
  /* dimensions of the inner grid*/
  
  int gridx2, gridy2, gridz2;   
  /* dimensions of the big grid*/
    
  int nr_energrads;
  /*contains all EnerGrads (either std or torque), for both biggrid and innergrid */  
  int shm_energrads; /*shared memory segment key: -1 by default*/
  
  //Defined is either this: 
  EnerGradStd *energrads_std; 
  //or this:
  EnerGradTorque *energrads_torque; 
  //or neither (if we calculate no potentials)
  
  int nr_neighbours;
  /* contains all Neighbours of the innergrid */
  int shm_neighbours; /*shared memory segment key: -1 by default*/
  Neighbour *neighbours;
    
  Potential *biggrid; 
  /* The outer grid (double grid spacing)
    It extends (gridextension*gridspacing = 32*0.9 A = 29 A)
     29 A beyond the inner grid)  
    The inner biggrid voxels are empty (NULL)
    
    Grid.biggrid is array of pointers-to-EnerGrid
    The pointers are internal: 
     they always point to a location in Grid.potentials*/   
  
  Voxel *innergrid;    
  /* The inner grid (single gridspacing)
     It covers a box around the protein and distcutoff (10.8 A) beyond it
       
     A Voxel contains a Potential (list of pointers-to-EnerGrid), 
     pointer-to-Neighbour and number of Neighbours
    
     Again, all of these pointers are internal:
       they point to locations in Grid.potentials and Grid.neighbours           
  */			    

  void init(double gridspacing0, int gridextension0, 
   double plateaudis0,double neighbourdis0,bool (&alphabet0)[MAXATOMTYPES]);
  void read_std(const char *filename);
  void read_torque(const char *filename);
  void write_std(const char *filename);
  void write_torque(const char *filename);

  double *_ratio; /*precomputed ratios to scale down distances, to reach 5 A*/  
  inline double get_ratio(double dissq) const {
    return _ratio[int(dissq*10000)];
  }  
  
  
  void calculate_std(int cartstatehandle, double plateaudis, double neighbourdis, float gridspacing, int gridextension, bool (&alphabet)[MAXATOMTYPES], bool cdie, float epsilon, bool calc_pot);
  void calculate_torque(int cartstatehandle, double plateaudis, double neighbourdis, float gridspacing, int gridextension, bool (&alphabet)[MAXATOMTYPES], bool cdie, float epsilon, bool calc_pot);
};

#endif
