#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "grid.h"

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

extern bool exists(const char *f);

extern void get_shm_name(int shm_id, char *shm_name);

int new_shm_id() {
  int ret = 1;
  char shm_name[100];
  while (1) {
    get_shm_name(ret, shm_name);
    int result = shm_open(shm_name, (O_CREAT|O_EXCL), (S_IREAD | S_IWRITE));
    if (result != -1) {
      return ret;
    }
    ret++;
  }
}

int main(int argc, char*argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Wrong number of arguments\n");
    fprintf(stderr, "Please provide grid file (input), grid header file (output)\n");
    return -1;
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "Grid file %s does not exist\n", argv[1]);
    fprintf(stderr, "Please provide grid file (input), grid header file (output)\n");
    return -1;    
  }

  Grid g;
  bool torquegrid;
#ifdef TORQUEGRID
  #define READ g.read_torque
  #define WRITE g.write_torque
#else
  #define READ g.read_std
  #define WRITE g.write_std
#endif  
  READ(argv[1]); 
  g.shm_energrads = new_shm_id();
  g.shm_neighbours = new_shm_id();
  printf("#shm %d %d\n", g.shm_energrads, g.shm_neighbours);
  WRITE(argv[2]);
}  
