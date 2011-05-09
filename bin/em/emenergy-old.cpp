#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <memory.h>

#ifdef GRADCHECK
const double d = 0.00001;
const double delta[][3] = {{d,0,0},{-d,0,0},{0,d,0},{0,-d,0},{0,0,d},{0,0,-d}};
#endif

#ifdef GRADCHECK2
const double d = 0.0001;
const double delta[][3] = {{d,0,0},{-d,0,0},{0,d,0},{0,-d,0},{0,0,d},{0,0,-d}};
#endif


double sigma_threshold = 3; //only consider voxels within 3 sigma


struct Map {
  double emweight;      // weight of the energy of this map 
  double overlapmargin; //fraction of maximum overlap that must be reached;less will give an energy penalty

  double situs_origx, situs_origy, situs_origz, situs_width;
  float resolution;
  double sigma;
  double base;
  double fac;

  float voxels_per_sigma;
  double maxdist;
  double maxdistsq;
  unsigned int dimx;
  unsigned int dimy;
  unsigned int dimz;
  double *densities;
  double *maxatoms;
  int nrmaxatoms;
  double maxsaturation;
  double maxsatweight;
  double *saturations;
  double *gradsatx;
  double *gradsaty;
  double *gradsatz;  
};

Map maps[100];
int nrmaps = 0;

typedef void (*Calcfunc) (Map &m, double dissq, double dx, double dy, double dz, double &totoverlap, int index, double &gradx, double &grady, double &gradz, bool &sat);

void calc_energy(Map &m, double dissq, double dx, double dy, double dz, double &totoverlap, int index, double &gradx, double &grady, double &gradz, bool &sat) {
  double density = m.densities[index];
  double overlap = density*pow(m.base, -dissq);
  totoverlap += overlap;
  /*
  double dis = sqrt(dissq);
  double grad = fac * -2 * dis * overlap;
  double fgrad = grad/dis;
  */
  double fgrad = m.fac * -2 * overlap;
  if (sat) {
    m.saturations[index] += overlap;
    double gx = dx*fgrad;
    double gy = dy*fgrad;
    double gz = dz*fgrad;
    gradx += gx;
    m.gradsatx[index] += gx;
    grady += gy;
    m.gradsaty[index] += gy;    
    gradz += gz;
    m.gradsatz[index] += gz;    
  }
  else {
    gradx += dx*fgrad;
    grady += dy*fgrad;
    gradz += dz*fgrad; 
  }
}

/*
void calc_overlap(Map &m, double dissq, double dx, double dy, double dz, double &totoverlap, int index, double &gradx, double &grady, double &gradz, bool &sat) {
  double oversaturation = m.saturations[index];
  if (oversaturation <= 0) return;
  double density = m.densities[index];
  //double saturation = oversaturation + density * m.maxsaturation;
  //double slope = -2 * m.maxsatweight * oversaturation/density;
  
  double slope = 1;
     
  double overlap = density*pow(m.base, -dissq);
  double fgrad = m.fac * -2 * overlap;
  double gradsat = m.gradsat
  gradx += dx*fgrad  * m.gradsatx[index]*density * slope;
  grady += dy*fgrad  * m.gradsaty[index]*density * slope;
  gradz += dz*fgrad  * m.gradsatz[index]*density * slope;
  
}
*/

unsigned int nratoms; 

extern "C" void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi);


inline void callfunc(
  Calcfunc func,
  
  int x, int y, int z, 
  double &ax, double &ay, double &az,

  Map &m, double dissq, double dx, double dy, double dz, double &totoverlap, int   index, double &gradx, double &grady, double &gradz, bool &sat
  ) 
  {

  #ifdef DISTCHECK	  
  double dissq_test = (x-ax)*(x-ax)+(y-ay)*(y-ay)+(z-az)*(z-az);   
  if (fabs(dissq-dissq_test) > 0.00001) 
  printf("%.5f %.5f %.5f\n",dissq-dissq_test, dissq, dissq_test);
  #endif
  
  #ifndef GRADCHECK2
  func(m, dissq, dx,dy,dz,totoverlap, index,gradx,grady,gradz, sat);  
  #else
  double overlap, overlap_new, grad;
  
  double totoverlap_old = totoverlap;
  double dgradx=0,dgrady=0,dgradz=0;
  func(m,dissq,dx,dy,dz,totoverlap,index,dgradx,dgrady,dgradz,sat);
  
  grad = sqrt(dgradx*dgradx+dgrady*dgrady+dgradz*dgradz);
  overlap =totoverlap_old - totoverlap;
  
  gradx += dgradx; grady+=dgrady;gradz+=dgradz;
  
  
  if (overlap == 0) return;
  
  double totoverlap0, dgradx0=0,dgrady0=0,dgradz0=0;
  
  double dissqnew;
  dissqnew = sqrt(dissq)+d;
  dissqnew *= dissqnew;
  totoverlap0 = totoverlap_old;
  bool fals = 0;
  func(m,dissqnew,dx,dy,dz,totoverlap0,index,dgradx0,dgrady0,dgradz0,fals);

  overlap_new = totoverlap_old - totoverlap0;
  
  printf("D+ %.3g %.3g %.3g %.3g %.3g\n", overlap, overlap_new, (overlap_new - overlap)/d, grad,  grad / ((overlap_new - overlap)/d));

  dissqnew = sqrt(dissq)-d;
  dissqnew *= dissqnew;
  totoverlap0 = totoverlap_old;
  func(m,dissqnew,dx,dy,dz,totoverlap0,index,dgradx0,dgrady0,dgradz0,fals);

  overlap_new = totoverlap_old - totoverlap0;
  
  printf("D- %.3g %.3g %.3g %.3g %.3g\n", overlap, overlap_new, (overlap_new - overlap)/-d, grad,  grad / ((overlap_new - overlap)/-d));
  
  
  for (int n=0;n<6;n++) {
    double dx2 = dx + delta[n][0];
    double dy2 = dy + delta[n][1];
    double dz2 = dz + delta[n][2];
    double dissq2 = dx2*dx2+dy2*dy2+dz2*dz2;
    totoverlap0 = totoverlap_old;
    func(m,dissq2,dx2,dy2,dz2,totoverlap0,index,dgradx0,dgrady0,dgradz0,fals);
    overlap_new = totoverlap_old - totoverlap0;    
    double dgrad = 0;
    if (delta[n][0] != 0) dgrad = dgradx;
    if (delta[n][1] != 0) dgrad = dgrady;
    if (delta[n][2] != 0) dgrad = dgradz;
    printf("delta%d %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n", n+1, dissq, dissq2, overlap, overlap_new, (overlap_new - overlap)/d, dgrad,  dgrad / ((overlap_new - overlap)/d));      
  }
  #endif
}


void iterate_voxels(Calcfunc func, Map &m, double ax, double ay, double az, 
double &totoverlap, double &gradx, double &grady, double &gradz) 
{
  bool sat = 0;
  if (m.maxsaturation > 0) sat = 1;
  
  int px0 = floor(ax);
  if (px0 < 0) px0=0;  
  if (px0 >= m.dimx) px0 = m.dimx-1;
  double dx0 = ax - px0;
  double dx0sq = dx0*dx0;
  int px1 = ceil(ax);
  if (px1 < 0) px1=0;  
  if (px1 == px0) px1++;
  if (px1 >= m.dimx) px1 = m.dimx-1;
  double dx1 = px1-ax;
  double dx1sq = dx1*dx1;

  int py0 = floor(ay);
  if (py0 < 0) py0=0;  
  if (py0 >= m.dimy) py0 = m.dimy-1;
  double dy0 = ay - py0;
  double dy0sq = dy0*dy0;
  int py1 = ceil(ay);
  if (py1 < 0) py1=0;    
  if (py1 == py0) py1++;    
  if (py1 >= m.dimy) py1 = m.dimy-1;
  double dy1 = py1-ay;
  double dy1sq = dy1*dy1;

  int pz0 = floor(az);
  if (pz0 < 0) pz0=0;  
  if (pz0 >= m.dimz) pz0 = m.dimz-1;
  double dz0 = az - pz0;
  double dz0sq = dz0*dz0;
  int pz1 = ceil(az);
  if (pz1 < 0) pz1=0;     
  if (pz1 == pz0) pz1++;    
  if (pz1 >= m.dimz) pz1 = m.dimz-1;
  double dz1 = pz1-az;
  double dz1sq = dz1*dz1;

  int indx,indxy, indxyz;
  const int dimxy = m.dimx*m.dimy;
  const int &dimx = m.dimx;

  //positive x loop
  indx = px1;
  double dsqx = dx1sq;
  double dx = dx1;
  int x;
  for (x = px1; x < m.dimx; x++) {
    if (dsqx > m.maxdistsq) break;

    //positive y loop
    indxy = indx+dimx*py1;
    double dsqxy = dsqx + dy1sq;
    double dy = dy1;
    int y;
    for (y = py1; y < m.dimy; y++){
      if (dsqxy > m.maxdistsq) break;

      //positive z loop
      indxyz = indxy+dimxy*pz1;
      double dsqxyz = dsqxy + dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < m.dimz; z++){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 
	  m, dsqxyz, dx,dy,dz,totoverlap, indxyz,gradx,grady,gradz, sat);

	//increment positive z loop
	indxyz += dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //negative z loop
      indxyz = indxy+dimxy*pz0;      
      dsqxyz = dsqxy + dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, dx,dy,-dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment negative z loop
	indxyz -= dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //increment positive y loop
      indxy += dimx;
      dsqxy += 2 *dy + 1;
      dy++;	
    }


    //negative y loop
    indxy = indx+dimx*py0;    
    dsqxy = dsqx + dy0sq;
    dy = dy0;
    for (y = py0; y >= 0; y--){
      if (dsqxy > m.maxdistsq) break;

      //positive z loop
      indxyz = indxy+dimxy*pz1;
      double dsqxyz = dsqxy + dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < m.dimz; z++){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, dx,-dy,dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment positive z loop
	indxyz += dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //negative z loop
      indxyz = indxy+dimxy*pz0;
      dsqxyz = dsqxy + dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, dx,-dy,-dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment negative z loop
	indxyz -= dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //increment negative y loop
      indxy -= dimx;
      dsqxy += 2 *dy + 1;
      dy++;	
    }


    //increment positive x loop
    indx++;
    dsqx += 2 *dx + 1;
    dx++;
  }

  //negative x loop
  indx = px0;  
  dsqx = dx0sq;
  dx = dx0;
  for (x = px0; x >= 0; x--) {
    if (dsqx > m.maxdistsq) break;

    //positive y loop
    indxy = indx+dimx*py1;        
    double dsqxy = dsqx + dy1sq;
    double dy = dy1;
    int y;
    for (y = py1; y < m.dimy; y++){
      if (dsqxy > m.maxdistsq) break;

      //positive z loop
      indxyz = indxy+dimxy*pz1;
      double dsqxyz = dsqxy + dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < m.dimz; z++){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, -dx,dy,dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment positive z loop
	indxyz += dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //negative z loop
      indxyz = indxy+dimxy*pz0;
      dsqxyz = dsqxy + dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, -dx,dy,-dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment negative z loop
	indxyz -= dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //increment positive y loop
      indxy += dimx;
      dsqxy += 2 *dy + 1;
      dy++;	
    }


    //negative y loop
    indxy = indx+dimx*py0;    
    dsqxy = dsqx + dy0sq;
    dy = dy0;
    for (y = py0; y >= 0; y--){
      if (dsqxy > m.maxdistsq) break;

      //positive z loop
      indxyz = indxy+dimxy*pz1;
      double dsqxyz = dsqxy + dz1sq;
      double dz = dz1;
      int z;
      for (z = pz1; z < m.dimz; z++){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, -dx,-dy,dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment positive z loop
	indxyz += dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //negative z loop
      indxyz = indxy+dimxy*pz0;
      dsqxyz = dsqxy + dz0sq;
      dz = dz0;
      for (z = pz0; z >= 0; z--){
        if (dsqxyz > m.maxdistsq) break;
	//dsqxyz is the actual distance squared

	callfunc(func,x,y,z,ax,ay,az, 	
	  m, dsqxyz, -dx,-dy,-dz,totoverlap, indxyz,gradx,grady,gradz,sat);

	//increment negative z loop
	indxyz -= dimxy;
	dsqxyz += 2 * dz + 1;
	dz++;
      }

      //increment negative y loop
      indxy -= dimx;
      dsqxy += 2 *dy + 1;
      dy++;	
    }


    //increment negative x loop
    indx--;    
    dsqx += 2 *dx + 1;
    dx++;
  }

}


/*
double situs_origx, situs_origy, situs_origz, situs_width;
*/

extern "C" void read_densitymaps_(char *densitymapsfile0, int len_densitymapsfile) {
  char buf[2000];
  
  char *densitymapsfile = new char[len_densitymapsfile];
  memcpy(densitymapsfile, densitymapsfile0, len_densitymapsfile);
  
  FILE *fil = fopen(densitymapsfile, "r");
  if (fil == NULL) {
    fprintf(stderr, "ERROR: file %s does not exist\n", densitymapsfile);
    exit(1);
  }
  
  while (!feof(fil)) {
    char densitymap[2000];
    char densitymax[2000];
    float resolution_value;
    if (!fgets(buf, 2000, fil)) break;
    if (buf[0] == '#') {
      buf[0] = ' ';
      continue;
    }
    float overlapmargin, emweight, maxsaturation,maxsatweight;
    int read = sscanf(buf, "%f %s %s %f %f %f %f", &resolution_value, densitymap, densitymax, &overlapmargin, &emweight, &maxsaturation, &maxsatweight);
    if (read <= 1) break;
    if (read != 7) {
      fprintf(stderr, "Reading error in %s:\n%s\n", densitymapsfile, buf);
      exit(1);
    }
    Map &m = maps[nrmaps];
    m.overlapmargin = overlapmargin;
    m.emweight = emweight;
    m.maxsaturation = maxsaturation;
    m.maxsatweight = maxsatweight;
    m.resolution = resolution_value;
    m.sigma = 0.5 * m.resolution / sqrt(3);
    double sigmafactor = 1/(2*m.sigma*m.sigma);
    m.fac = (sigmafactor*sigmafactor)/(sigmafactor+sigmafactor);
    m.base = exp(m.fac);

    read_vol(densitymap, &m.situs_width, &m.situs_origx,  &m.situs_origy,&m.situs_origz,&m.dimx,&m.dimy,&m.dimz,&m.densities);

    const int maxmaxatoms = 100000;
    double maxatoms[maxmaxatoms];
    FILE *fil2 = fopen(densitymax, "r");
    int nrmaxatoms = 0;
    while(!feof(fil2)) {
      if (!fgets(buf, 2000, fil2)) break;
      float max;
      if (!sscanf(buf, "%*s %*d %*s %f", &max)) break;
      if (nrmaxatoms == maxmaxatoms) {
        fprintf(stderr, "%s contains more than %d atoms\n", densitymax, maxmaxatoms);
	exit(1);
      }
      maxatoms[nrmaxatoms] = max;
      nrmaxatoms++;      
    }
    m.maxatoms = new double[nrmaxatoms];
    memcpy(m.maxatoms, maxatoms, nrmaxatoms*sizeof(double));
    m.nrmaxatoms = nrmaxatoms;
    m.voxels_per_sigma = m.sigma / m.situs_width;
    m.maxdist = sigma_threshold * m.voxels_per_sigma;
    m.maxdistsq = m.maxdist * m.maxdist;    
    
    if (m.maxsaturation > 0) {    
     int gridsize = m.dimx * m.dimy * m.dimz;      
      m.saturations = new double[gridsize];
      memset(m.saturations, 0,gridsize*sizeof(double));
      m.gradsatx = new double[gridsize];
      memset(m.gradsatx, 0,gridsize*sizeof(double));
      m.gradsaty = new double[gridsize];
      memset(m.gradsaty, 0,gridsize*sizeof(double));
      m.gradsatz = new double[gridsize];
      memset(m.gradsatz, 0,gridsize*sizeof(double));    
    }    
    nrmaps++;
  }
  delete[] densitymapsfile;
}

double emenergy (Map &m, int &nratoms, double *atoms, double *forces, int &update_forces, double &totgradx, double &totgrady, double &totgradz) {  
  double energy = 0;
  if (nratoms != m.nrmaxatoms) {
    int mapnr = &m - maps+1;
    fprintf(stderr, "ERROR: Electron density map %d: number of atoms in PDB file (%d) is not equal in the number of atoms in the maximum density file (%d)\n", mapnr, nratoms, m.nrmaxatoms);
    exit(0);
  }
  if (m.maxsaturation > 0) {
    int gridsize = m.dimx * m.dimy * m.dimz;      
    memset(m.saturations, 0,gridsize*sizeof(double));
    memset(m.gradsatx, 0,gridsize*sizeof(double));
    memset(m.gradsaty, 0,gridsize*sizeof(double));
    memset(m.gradsatz, 0,gridsize*sizeof(double));    
  }
  totgradx = 0; totgrady = 0; totgradz = 0;
    
  for (int n=0;n<nratoms;n++) {
    double ax = (atoms[6*n] - m.situs_origx)/m.situs_width;    
    double ay = (atoms[6*n+2] - m.situs_origy)/m.situs_width;
    double az = (atoms[6*n+4] - m.situs_origz)/m.situs_width;
    double totoverlap = 0;
    double gradx = 0, grady = 0, gradz = 0;

    iterate_voxels(calc_energy, m, ax,ay,az,totoverlap,gradx,grady,gradz);
  
    double max_overlap = m.maxatoms[n];  
    double gradscale = 2*m.emweight/(-max_overlap*m.situs_width);
    totoverlap /= max_overlap;
    double curr_energy = 0;
    if (totoverlap < m.overlapmargin) {
      double lack = m.overlapmargin - totoverlap;
      curr_energy = m.emweight * lack*lack;
      energy += curr_energy;
      double currgradscale = lack * gradscale;
      if (update_forces == 1) {
        /*
        gradx /= -max_overlap*m.situs_width;
        gradx *= 2 * emweight * lack;
        grady /= -max_overlap*m.situs_width;
        grady *= 2 * emweight * lack;
        gradz /= -max_overlap*m.situs_width; 
        gradz *= 2 * emweight * lack;
        */
        gradx *= currgradscale;
        grady *= currgradscale;
        gradz *= currgradscale;
	totgradx += gradx; totgrady += grady; totgradz += gradz;	
        forces[3*n] += gradx;
        forces[3*n+1] += grady;
        forces[3*n+2] += gradz;
      }
    }
  }
  if (m.maxsaturation > 0) {
    double clashenergy = 0;
    int gridsize = m.dimx * m.dimy * m.dimz;      
    for (int n=0; n < gridsize; n++) {
      double d = m.densities[n];
      if (d == 0) continue;
      double sat = m.saturations[n]/d;
      double excess = sat - m.maxsaturation;
      if (excess <= 0) continue;
      excess /= m.maxsaturation;
      clashenergy += (excess*excess)*m.maxsatweight;
      
      if (update_forces == 1) {
        double slope = -2*excess * m.maxsatweight;
	double gradx = m.gradsatx[n]/d * slope/m.maxsaturation;
	double grady = m.gradsaty[n]/d * slope/m.maxsaturation;
	double gradz = m.gradsatz[n]/d * slope/m.maxsaturation;	
        gradx /= m.situs_width;
        grady /= m.situs_width;
        gradz /= m.situs_width;	
	totgradx += gradx; totgrady += grady; totgradz += gradz;		
      }
      
    }
    printf("EM energy: %d %.3f %.3f\n", &m-maps+1, energy, clashenergy);  
    energy += clashenergy;
    if (update_forces) {
      double forcex = totgradx/nratoms;
      double forcey = totgrady/nratoms;
      double forcez = totgradz/nratoms;      
      for (int n=0;n<nratoms;n++) {
        /*      
        double ax = (atoms[6*n] - m.situs_origx)/m.situs_width;    
        double ay = (atoms[6*n+2] - m.situs_origy)/m.situs_width;
        double az = (atoms[6*n+4] - m.situs_origz)/m.situs_width;
        double gradx = 0, grady = 0, gradz = 0; double bogus;
        iterate_voxels(calc_overlap, m, ax,ay,az,bogus,gradx,grady,gradz); 
        gradx /= m.situs_width;
        grady /= m.situs_width;
        gradz /= m.situs_width;	
	totgradx += gradx; totgrady += grady; totgradz += gradz;	
        forces[3*n] += gradx;
        forces[3*n+1] += grady;
        forces[3*n+2] += gradz;
	*/
        forces[3*n] += forcex;
        forces[3*n+1] += forcey;
        forces[3*n+2] += forcez;	
      }
    } 
  }
  else {
    printf("EM energy: %d %.3f\n", &m-maps+1, energy);  
  }
  
  return energy;  
}


extern "C" void emenergy_(double &energy, int &nratoms, double *atoms, double *forces, int &update_forces) {  

  energy = 0;
  for (int m = 0; m < nrmaps; m++) {

#ifdef GRADCHECK
  
    double tgradx,tgrady,tgradz, tgradx0,tgrady0,tgradz0;
    double denergy = emenergy(maps[m], nratoms, atoms, forces, update_forces, tgradx,tgrady,tgradz);
    energy += denergy;
    
    double *datoms = new double[6*nratoms];
    memcpy(datoms, atoms, 6*nratoms*sizeof(double));
    double energy0;
    //printf("tgradx %.6g tgrady %.6g tgradz %.6g\n", tgradx,tgrady,tgradz);
    for (int i=0;i<6;i++) {
      memcpy(datoms, atoms, 6*nratoms*sizeof(double));
      for (int n=0;n<nratoms;n++) {
        datoms[6*n] += delta[i][0];
        datoms[6*n+2] += delta[i][1];
        datoms[6*n+4] += delta[i][2];
      }    
      double grad = 0;
      if (delta[i][0] != 0) grad = tgradx;
      if (delta[i][1] != 0) grad = tgrady;
      if (delta[i][2] != 0) grad = tgradz;
      int no_update = 0;
      energy0 = emenergy(maps[m], nratoms, datoms, forces, no_update, tgradx0,tgrady0,tgradz0);
      printf("trans%d %.3f %.6g %.6g %.6g\n", i, energy0, (energy0-denergy)/d, grad, grad/((energy0-denergy)/d));
    }
    printf("\n");
  }
}

#else
    //Normal execution
    
    double dumx,dumy,dumz;
    double denergy = emenergy(maps[m], nratoms, atoms, forces, update_forces,dumx,dumy,dumz);
    energy += denergy;
  }
  #ifdef GRADCHECK2
  exit(0);
  #endif
}

#endif
