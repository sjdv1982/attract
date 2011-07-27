#include <cstdio>
#include <cstring>
#include <sys/mman.h>

char *files[1000];
int nrfiles;

void get_shm_name(int shm_id, char *shm_name) {
  sprintf(shm_name, "/attract-grid%d", shm_id);
}

int main(int argc, char *argv[]) {  
  files[0] = "/attract-prox-100000-36000-200000-31";
  files[1] = "/attract-prox-100000-64000-200000-31";
  files[2] = "/attract-prox-100000-100000-200000-31";
  files[3] = "/attract-prox-100000-100000-200000-48";
  files[4] = "/attract-prox-100000-144000-200000-48";  
  //files[1] = "/attract-prox-250";
  for (int n = 1; n <= 100; n++) {
    char *name = new char[100];
    get_shm_name(n,name);
    files[n+1] = name;
  }
  nrfiles = 101;
  for (int n = 0; n < nrfiles; n++) {
    //printf("%s\n", files[n]);
    shm_unlink(files[n]);
  }
}
