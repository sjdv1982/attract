#include <cstdio>
#include <cstring>
#include <sys/mman.h>

char *files[1000];
int nrfiles;

void get_shm_name(int shm_id, char *shm_name) {
  sprintf(shm_name, "/attract-grid%d", shm_id);
}

int main(int argc, char *argv[]) {  
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
