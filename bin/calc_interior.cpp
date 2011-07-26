
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <cstring>


#include "makegrid.h"

//A grid is constructed 10.8 A around the protein
//Everything within 8.5 A of at least 91 atoms is assigned as interior, and holes are filled
//finally, the interior is shrunk by 9 voxels

typedef struct {
	char  recd[7];        /*       1 -  6 */
	int   serial;         /*       7 - 11 */
	char  type[3];        /*      13 - 14 */
	char  loc[3];         /*      15 - 16 */
	char  alt[2];         /*           17 */
	char  res[5];         /*      18 - 21 */
	char  chain[2];       /*           22 */
	int   seq;            /*      23 - 26 */
	char  icode[2];       /*           27 */
	float x;              /*      31 - 38 */
	float y;              /*      39 - 46 */
	float z;              /*      47 - 54 */
	float occupancy;      /*      55 - 60 */
	float beta;           /*      61 - 66 */
	int   footnote;       /*      68 - 70 */
	char  segid[5];       /*      73 - 76 */
	char  element[3];     /*      77 - 78 */
	char  charge[3];      /*      79 - 80 */
	float weight;         /* mass of atom */
} PDB;

extern "C" void read_pdb(char *file_name, unsigned *num_atoms, PDB **in_pdb);

extern "C" void write_vol(char *vol_file, double width, double origx, double origy, double origz, unsigned extx, unsigned exty, unsigned extz, double *phi);

using namespace std;

typedef double Coor[3];


inline void inc(Coor *dis, int nrdis, int index, double d) {
  for (int n = 0; n < nrdis; n++) 
    dis[n][index] += d;
}

inline bool test_surface(Coor *dis, int nrdis, int inside_threshold) {
  int count = 0;
  for (int n = 0; n < nrdis; n++) {
    Coor &d = dis[n];
    double dsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
    if (dsq < interior_cutoffsq) {
      count++;
      if (count >= inside_threshold) return 1;
    }
  }
  return 0;
}

int scan_count = 0;

inline bool scan_line(double *grid, int index0, int dim, int indincrement) {
  bool change = 0;
  int mode = 1; //0 = searching, 1 = filling
  int index = index0;
  for (int i = 0; i < dim; i++) {
    double &v = grid[index];
    index += indincrement;
    if (mode == 0) {
      if (v == 1) mode = 1;
    }
    else if (mode == 1) {
      if (v == 0) {
        v = 1;
	scan_count++;
	change = 1;
      } 
      else if (v == interior_value) {
        mode = 0;
      }
    }
    
  }
  mode = 1;
  index = index0 + (dim-1) * indincrement;
  for (int i = dim-1; i >= 0; i--) {
    double &v = grid[index];
    index -= indincrement;
    if (mode == 0) {
      if (v == 1) mode = 1;
    }
    else if (mode == 1) {
      if (v == 0) {
        v = 1;
	scan_count++;
	change = 1;
      } 
      else if (v == interior_value) {
        mode = 0;
      }
    }
    
  }
    
  return change;
}


void scan_interior(double *grid, unsigned int gridx, unsigned int gridy, unsigned int gridz) {
  bool change;
  int n,nn;
  int iteration = 0;
  do {
    change = 0;
    iteration++;
    //scan X lines
    for (n = 0; n < gridy; n++) {
      for (nn = 0; nn < gridz; nn++) {
        int index = gridx*n+gridx*gridy*nn;
        if (scan_line(grid, index, gridx, 1)) change = 1;
      }
    }

    //scan Y lines
    for (n = 0; n < gridx; n++) {
      for (nn = 0; nn < gridz; nn++) {
        int index = n+gridx*gridy*nn;
        if (scan_line(grid, index, gridy, gridx)) change = 1;
      }
    }

    //scan Z lines
    for (n = 0; n < gridx; n++) {
      for (nn = 0; nn < gridy; nn++) {
        int index = n+gridx*nn;
        if (scan_line(grid, index, gridz, gridx*gridy)) change = 1;
      }
    }
    printf("Scan %d %d\n", iteration, scan_count);

  } while (change);
}
int main(int argc, char*argv[]) {
  unsigned int n;
  if (argc < 2) {
    printf("Too few arguments\n");
    printf("Please PDB file, SITUS output file\n");
    return -1;
  }
  PDB *pdb;
  unsigned int nratoms;
  read_pdb(argv[1], &nratoms, &pdb);
   
  float minx=99999, miny=99999, minz=99999;
  float maxx=-99999, maxy=-99999, maxz=-99999;
  for (n = 0; n < nratoms; n++) {
    PDB &a = pdb[n];
    if (a.x > maxx) maxx = a.x; if (a.x < minx) minx = a.x;
    if (a.y > maxy) maxy = a.y; if (a.y < miny) miny = a.y;
    if (a.z > maxz) maxz = a.z; if (a.z < minz) minz = a.z;    
  }
  
  printf("%.3f %.3f\n", minx, maxx);
  printf("%.3f %.3f\n", miny, maxy);
  printf("%.3f %.3f\n", minz, maxz);  
  
  unsigned int gridx = ceil((maxx-minx+2*boxspace)/gridspacing);
  unsigned int gridy = ceil((maxy-miny+2*boxspace)/gridspacing);
  unsigned int gridz = ceil((maxz-minz+2*boxspace)/gridspacing);
  
  printf("%d %d %d\n", gridx, gridy, gridz);
  
  double *grid = new double[gridx * gridy * gridz];
  memset(grid,0,long(gridx) * long(gridy) * long(gridz)*sizeof(double));
  double orix = minx-boxspace;
  double oriy = miny-boxspace;
  double oriz = minz-boxspace;
  
  Coor *dis = new Coor[nratoms];
  for (n = 0; n < nratoms; n++) {
    dis[n][0] = -(pdb[n].x-orix);
    dis[n][1] = -(pdb[n].y-oriy);
    dis[n][2] = -(pdb[n].z-oriz);
  }
  Coor *disx = dis;
  Coor *disxy = new Coor[nratoms];
  Coor *disxyz = new Coor[nratoms];
  
  int indx = 0, indxy = 0, indxyz = 0;
  int insidec = 0;
  for (int x = 0; x < gridx; x++) {        
    indxy = indx;
    memcpy(disxy, disx, nratoms*sizeof(Coor));
    for (int y = 0; y < gridy; y++) {
      
      indxyz = indxy;
      memcpy(disxyz, disxy, nratoms*sizeof(Coor));    
      for (int z = 0; z < gridz; z++) {         
	 bool t = test_surface(disxyz, nratoms, inside_threshold);
	 if (t) {
	   grid[indxyz] = interior_value;
	   insidec++;
	 }
	 inc(disxyz, nratoms,  2, gridspacing);
	 indxyz += gridx * gridy;
      }
      
      inc(disxy, nratoms,1, gridspacing);
      indxy += gridx;
      
    }
    
    inc(disx, nratoms,0, gridspacing);
    indx += 1;
    
  }      

  printf("inside: %d (%.3f %%)\n", insidec, float(long(100) * long(insidec))/(long(gridx)*long(gridy)*long(gridz)));	      
  scan_interior(grid, gridx, gridy, gridz);
  
  int outside = scan_count;
  printf("outside: %d (%.3f %%)\n", outside, float(long(100) * long(outside))/(long(gridx)*long(gridy)*long(gridz)));	      
  
  indx = 0; 
  for (int x = 0; x < gridx; x++) {        
    indxy = indx;
    for (int y = 0; y < gridy; y++) {      
      indxyz = indxy;
      for (int z = 0; z < gridz; z++) {         
         if (grid[indxyz] == 0) {
	   grid[indxyz] = interior_value;
	   insidec++;
	 }
	 indxyz += gridx * gridy;
      }      
      indxy += gridx;      
    }    
    indx += 1;    
  }      
  
  printf("inside, after scan: %d (%.3f %%)\n", insidec, float(long(100) * long(insidec))/(long(gridx)*long(gridy)*long(gridz)));	      
  
  
  if (shrink_interior > 0) {
    double *grid0 = new double[long(gridx)*long(gridy)*long(gridz)];
    memcpy(grid0, grid, long(gridx)*long(gridy)*long(gridz)*sizeof(double)); 

    //Z scan
    for (int x = 0; x < gridx; x++) {
      for (int y = 0; y < gridy; y++) {
        bool mode = 0;
	for (int z = 0; z < gridz; z++) {
          int index = x+gridx*y+gridx*gridy*z;
	  bool inside = (fabs(grid0[index] - interior_value) <0.1);
	  if (mode == 0 && inside == 1) {
	    int zz1 = z + shrink_interior;
	    if (zz1 >= gridz) zz1 = gridz - 1;
	    for (int zz = z; zz < zz1; zz++) {
	      int index2 = x+gridx*y+gridx*gridy*zz;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;
	      grid[index2] = 0;
	      
	    }
	    mode = 1;
	    continue;
	  }
	  if (mode == 1 && inside == 0) {
	    int zz0 = z - shrink_interior;
	    if (zz0 < 0) zz0 = 0;
	    for (int zz = zz0; zz < z; zz++) {
	      int index2 = x+gridx*y+gridx*gridy*zz;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;
	      grid[index2] = 0;	      
	    }
	    mode = 0;
	  }	  
        }
      }
    }
    
    //Y scan
    for (int x = 0; x < gridx; x++) {
      for (int z = 0; z < gridz; z++) {
        bool mode = 0;
	for (int y = 0; y < gridy; y++) {
          int index = x+gridx*y+gridx*gridy*z;
	  bool inside = (fabs(grid0[index] - interior_value) <0.1);
	  if (mode == 0 && inside == 1) {
	    int yy1 = y + shrink_interior;
	    if (yy1 >= gridy) yy1 = gridy - 1;
	    for (int yy = y; yy < yy1; yy++) {
	      int index2 = x+gridx*yy+gridx*gridy*z;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;
	      grid[index2] = 0;
	      
	    }
	    mode = 1;
	    continue;
	  }
	  if (mode == 1 && inside == 0) {
	    int yy0 = y - shrink_interior;
	    if (yy0 < 0) yy0 = 0;
	    for (int yy = yy0; yy < y; yy++) {
	      int index2 = x+gridx*yy+gridx*gridy*z;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;	      
	      grid[index2] = 0;
	    }
	    mode = 0;
	  }	  
        }
      }
    }

    //X scan
    for (int y = 0; y < gridy; y++) {
      for (int z = 0; z < gridz; z++) {
        bool mode = 0;
	for (int x = 0; x < gridx; x++) {
          int index = x+gridx*y+gridx*gridy*z;
	  bool inside = (fabs(grid0[index] - interior_value) <0.1);
	  if (mode == 0 && inside == 1) {
	    int xx1 = x + shrink_interior;
	    if (xx1 >= gridx) xx1 = gridx - 1;
	    for (int xx = x; xx < xx1; xx++) {
	      int index2 = xx+gridx*y+gridx*gridy*z;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;
	      grid[index2] = 0;
	      
	    }
	    mode = 1;
	    continue;
	  }
	  if (mode == 1 && inside == 0) {
	    int xx0 = x - shrink_interior;
	    if (xx0 < 0) xx0 = 0;
	    for (int xx = xx0; xx < x; xx++) {
	      int index2 = xx+gridx*y+gridx*gridy*z;
	      if (grid0[index2] == 0) printf("ERR\n");
	      if (grid[index2] == interior_value) insidec--;
	      grid[index2] = 0;	      
	    }
	    mode = 0;
	  }	  
        }
      }
    }
    
  }
  printf("inside, after grow/shrink: %d (%.3f %%)\n", insidec,  float(long(100) * long(insidec))/(long(gridx)*long(gridy)*long(gridz)));	      
  write_vol(argv[2], gridspacing, orix, oriy, oriz, gridx, gridy, gridz, grid);

}
