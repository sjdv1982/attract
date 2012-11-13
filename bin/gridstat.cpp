#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include  "max.h"


bool exists(const char *f) {
  FILE *fil = fopen(f, "r");
  if (fil == NULL) return 0;
  else {
    fclose(fil);
    return 1;
  }
}

typedef double Coor[3]; /*Coor is double[3]*/
typedef float Coorf[3]; /*Coorf is float[3]*/

struct EnerGrad {
  float energy;
  Coorf grad;
#ifdef TORQUEGRID  
  float torques[9];
#endif  
};


double gridspacing; /* 0.9 A */
int gridextension; /* 32; always make this number even!*/
double plateaudis; /* 10 A */
double plateaudissq;
double plateaudissqinv;
double neighbourdis; 
double neighbourdissq;
double proxlim, proxmax;
int proxmaxtype;
int natoms;
int nhm;
double *modedofs;
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
/*contains all EnerGrads, for both biggrid and innergrid */  
int shm_energrads; /*shared memory segment key: -1 by default*/

int nr_neighbours;

struct Neighbour {
  char type; /*1 = under 5 A, 2 = 5-7 A*/
  unsigned int index; /*receptor atom index*/
};

void error(const char *filename) {
  fprintf(stderr, "Reading error in grid file %s\n", filename);
  exit(1);
}

void do_read(const char *filename) {
  int read;
  
  FILE *f = fopen(filename, "rb");  
  if (f == NULL) error(filename);

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
  read = fread(&nhm,sizeof(int),1,f);
  if (!read) error(filename);
  if (nhm > 0) {
    modedofs = new double[nhm];
    read=fread(modedofs,nhm*sizeof(double),1,f);
    if (!read) error(filename);
  }
  read=fread(&pivot,sizeof(Coor),1,f);
  if (!read) error(filename);
  
  read = fread(&nr_energrads, sizeof(nr_energrads),1,f);  
  if (!read) error(filename);
  read = fread(&shm_energrads, sizeof(shm_energrads),1,f);  
  if (!read) error(filename);

  if (nr_energrads) {
    if (shm_energrads == -1) {
      int shift = nr_energrads*sizeof(EnerGrad);
      fseek(f, shift, SEEK_CUR);
      if (!read) error(filename);
    }
    else {
      error(filename);
    }  
  }

  read = fread(&nr_neighbours, sizeof(nr_neighbours),1,f);  
  if (!read) error(filename);
  
  printf("Number of potentials: %d Number of neighbours: %d\n", nr_energrads, nr_neighbours);
}

int main(int argc, char*argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Wrong number of arguments\n");
    fprintf(stderr, "Please provide grid file (input)\n");
    return -1;
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "Grid file %s does not exist\n", argv[1]);
    fprintf(stderr, "Please provide grid file (input)\n");
    return -1;    
  }
  do_read(argv[1]);
}  
