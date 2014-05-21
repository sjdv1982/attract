#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <memory>
using namespace std;

/*
  cluster_struc.cpp
  author: Sjoerd J. de Vries
  Clusters solutions based on pairwise RMSD matrix
   RMSD matrix format:
   <#structure1(range: 1-MAX_STRUC)> <#structure2> <RMSD>
    order must be:
    1 2, 1 3, .., 1 #struc, 2 1, 2 2, ..., #nstruc-1 #nstruc
  Structural neighbors are defined using a supplied cutoff
  Each time, the structure having the largest number of neighbors
   is taken as seed for a new cluster
  Only clusters larger than minsize are selected
  if the option -f is set, full linkage is applied
*/

#define MAX_STRUC 100000

void full_link(bool **neighbor , int currstrucnr, int nrstruc, bool *linked_neighbor, bool *done) {
  /*finds all neigbors of neighbors (and their neighbors etc.)
     and includes them in the neigbor list
    linked_neighbor contains the linked neigbor list
    done contains a list of structures that have already been
    examinated
  */
  static int recursivity;
  recursivity++;
  //printf("%d %d\n",currstrucnr,recursivity);
  done[currstrucnr] = 1;
  for (int n = 0; n < nrstruc; n++) {
    if (neighbor[currstrucnr][n] && !done[n]) {
      linked_neighbor[n] = 1;
      full_link(neighbor, n, nrstruc, linked_neighbor, done);
    }
  }
  recursivity--;
  return;
}

int CheckOption(char *arg, bool &full_linkage) {
  if (!strcmp(arg, "-f")) {
    full_linkage = 1;
  }
  if (arg[0] == '-') return 1;
  return 0;
}

int main(int argc, char *argv[]) {
  int n;
  float cutoff;
  int minsize;
  bool full_linkage = 0;

 
//Check command line
  if (argc < 4) {
    fprintf(stderr, "Wrong number of arguments\n");
    fprintf(stderr, "Usage: cluster_struc [-f] <RMSD file> <cutoff> <minsize>\n");
    return 1;
  }

  int optioncounter = 0;
  bool res;

  res = 1;
  do {
    res = CheckOption(argv[optioncounter+1], full_linkage);
    optioncounter += res;
  } while(res);
  int filearg = optioncounter + 1;

  FILE *f = fopen(argv[filearg], "r");
  if (!f) {
    fprintf(stderr, "File %s does not exist\n", argv[filearg]);
    fprintf(stderr, "Usage: cluster_struc [-f] <RMSD file> <cutoff> <minsize>\n");
    return 2;
  }

  res = 1;
  do {
    res = CheckOption(argv[optioncounter+2], full_linkage);
    optioncounter += res;
  } while(res);

  cutoff = atof(argv[optioncounter + 2]);
  if (cutoff <= 0) {
    fprintf(stderr, "Neighbor cutoff %f out of range\n", argv[optioncounter + 2]);
    fprintf(stderr, "Usage: cluster_struc [-f] <RMSD file> <cutoff> <minsize>\n");
    return 3;
  }

  res = 1;
  do {
    res = CheckOption(argv[optioncounter+3], full_linkage);
    optioncounter += res;
  } while(res);


  minsize = atoi(argv[optioncounter + 3]);
  if (minsize < 0) {
    fprintf(stderr, "Minimal cluster size %d out of range\n", minsize);
    fprintf(stderr, "Usage: cluster_struc [-f] <RMSD file> <cutoff> <minsize>\n");
    return 4;
  }

  for (n = optioncounter + 4; n < argc; n++) {
    CheckOption(argv[n], full_linkage);
  }

//Determine number of structures
  char buf[1000];
  int nrstruc = 0;
  while (!feof(f)) {
    fgets(buf, 1000, f);
    int first, second;
    sscanf(buf, "%d %d", &first, &second);
    if (first != 1) {
      if (!nrstruc) {
        fprintf(stderr, "Reading error in %s\n", argv[filearg]);
        return 5;
      }
      nrstruc++;
      break;
    }
    nrstruc++;
    if (nrstruc >= MAX_STRUC) {
      fprintf(stderr, "Matrix %s is too large\n", argv[filearg]);
      return 6;
    }
  }

//Define boolean neighbor matrix
  int *neighborcount = new int[nrstruc];
  memset(neighborcount, 0, nrstruc * sizeof(int));
  bool **neighbor = new bool*[nrstruc];
  for (n = 0; n < nrstruc; n++) {
    neighbor[n] = new bool[nrstruc];
    memset(neighbor[n], 0, nrstruc * sizeof(bool));
  }

//Parse the file
  fclose(f);
  f = fopen(argv[filearg], "r");
  for (n = 0; n < nrstruc; n++) {
    for (int nn = n+1; nn < nrstruc; nn++) {
      bool err = 0;
      int first, second;
      float rmsd;
      do {
        if (!fgets(buf, 1000, f)) {
	  err = 1;
	  break;
	}
	if (sscanf(buf, "%d %d %f", &first, &second, &rmsd) != 3) {            err = 1;
	  break;
	}
        if (first != n+1 || second != nn+1) {
	  err = 1;
	  break;
	}
      } while (0);
      if (err) {
        fprintf(stderr, "Reading error in %s\n", argv[filearg]);
        for (int i = 0; i < nrstruc; i++) {
          delete[] neighbor[i];
        }
        delete[] neighbor, neighborcount;
        return 5;
      }
      if (rmsd < cutoff) {
        neighbor[n][nn] = 1;
        neighbor[nn][n] = 1;
        neighborcount[n]++;
	neighborcount[nn]++;
      }
    }

  }

//Clustering cycle
  int clustercounter = 0;
  do {
  //Select the largest cluster
    int cluster_max = -1;
    int cluster_nr = -1;
    for (n = 0; n < nrstruc; n++) {
      if (neighborcount[n] > cluster_max) {
        cluster_max = neighborcount[n];
	cluster_nr = n;
      }
    }

  //Terminate if no suitably large clusters are found
    if (cluster_max < minsize -1) break;

  //Print the cluster
    bool *clust;bool *done;
    if (full_linkage) {
      clust = new bool[nrstruc];
      memset(clust, 0, nrstruc * sizeof(bool));
      done = new bool[nrstruc];
      memset(done, 0, nrstruc * sizeof(bool));
      full_link(neighbor, cluster_nr, nrstruc, clust, done);
    }
    else {
      clust = neighbor[cluster_nr];
    }

    clustercounter++;
    printf("Cluster %d -> %d", clustercounter, cluster_nr + 1);
    for (n = 0; n < nrstruc; n++) {
      if (clust[n]) printf(" %d", n+1);
    }
    printf("\n");

  //Eliminate its members
    neighborcount[cluster_nr] = -1;
    for (n = 0; n < nrstruc; n++) {
      if (clust[n]) {
        neighborcount[n] = -1;
        for (int nn = 0; nn < nrstruc; nn++) {
	  if (neighbor[nn][n]) neighborcount[nn]--;
	  neighbor[nn][n] = 0;
	}
      }
    }
    if (full_linkage) {
      delete [] clust;
      delete [] done;
    }

  } while (1);

//Free memory
  for (n = 0; n < nrstruc; n++) {
    delete[] neighbor[n];
  }
  delete[] neighbor;

  return 0;
}
