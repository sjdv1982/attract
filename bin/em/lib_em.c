/* Copyright 2012 Sjoerd de Vries
this code is licensed under the GPL version 2 or later
contains code from pdb2vol from Situs, copyright Willy Wriggers
*/

#include "lib_em.h"

double gridify (
 Coor *pdb2, int num2, 
 double *phi2,  int nvox2, 
 double width, int g_extx, int g_exty, int extz,
 double minx, double miny, double minz
) {
  double varp, mass_total;  
  int count, i;
  double gx,gy,gz,a,b,c, w, ai,bi,ci;
  int x0,y0,z0,x1,y1,z1;
  for (count=0;count<nvox2;count++) *(phi2+count) = 0.0;
  varp = 0.0;
  mass_total = 0.0;
  for (i=0;i<num2;++i) {
    /* compute position within grid*/
    gx = (pdb2[i][0] - minx) / width;
    gy = (pdb2[i][1] - miny) / width;
    gz = (pdb2[i][2] - minz) / width;
    x0 = floor (gx); 
    y0 = floor (gy); 
    z0 = floor (gz); 
    if (x0 >= g_extx) continue;
    if (y0 >= g_exty) continue;
    if (z0 >= extz) continue;
    x1 = x0+1;
    y1 = y0+1; 
    z1 = z0+1; 
    if (x1 < 0) continue;
    if (y1 < 0) continue;
    if (z1 < 0) continue;
    /* interpolate */
    a = x1-gx;
    if (x0 < 0) a = 0;
    b = y1-gy;
    if (y0 < 0) b = 0;
    c = z1-gz;
    if (z0 < 0) c = 0;
    ai = 1-a;
    if (x1 == g_extx) ai = 0;
    bi = 1-b;
    if (y1 == g_exty) bi = 0;        
    ci = 1-c;
    if (z1 == extz) ci = 0;        
    
    w = a * b * c;
    if (w > 0) { 
      *(phi2+g1idz(z0,y0,x0)) += w; varp += a * b * c * (ai*ai+bi*bi+ci*ci);
    }
    w = a * b * ci;
    if (w > 0) {
      *(phi2+g1idz(z1,y0,x0)) += w; varp += a * b * ci * (ai*ai+bi*bi+c*c);
    }
    w = a * bi * c;
    if (w > 0) {
      *(phi2+g1idz(z0,y1,x0)) += w; varp += a * bi * c * (ai*ai+b*b+ci*ci);
    }
    w = ai * b * c;
    if (w > 0) {
      *(phi2+g1idz(z0,y0,x1)) += w; varp += ai * b * c * (a*a+bi*bi+ci*ci);
    }      
    w = a * bi * ci;
    if (w > 0) {
      *(phi2+g1idz(z1,y1,x0)) += w; varp += a * bi * ci * (ai*ai+b*b+c*c);
    }      
    w = ai * bi * c;
    if (w > 0) {
      *(phi2+g1idz(z0,y1,x1)) += w; varp += ai * bi * c * (a*a+b*b+ci*ci);
    }      
    w = ai * b * ci;
    if (w > 0) {
      *(phi2+g1idz(z1,y0,x1)) += w; varp += ai * b * ci * (a*a+bi*bi+c*c);
    }      
    w = ai * bi * ci;
    if (w > 0) {
      *(phi2+g1idz(z1,y1,x1)) += w; varp += ai * bi * ci * (a*a+b*b+c*c);
    }	
    mass_total ++;
  }
  varp /= mass_total;
  return varp;
}

double *calc_kernel(double width, double reso, double kampl, unsigned int *pg_extx) {
    double *phi;
    int exth, g_extx, g_exty, extz, nvox;
    double bvalue, cvalue;
    double rh1, rc1, sig1, sigmap, dsqu, kmsd, varmap;
    int count, indx, indy, indz;
       
    sig1 = reso/2.0;
    rh1 = sig1 * sqrt(log(2.0)) / sqrt(1.5);
    rc1 = sqrt(3.0)*sig1;
    kmsd = sig1*sig1/(width*width);
    varmap = kmsd; 
    
    sigmap = sqrt(varmap/3.0); /* sigma-1D */
    exth = (int) ceil(3*sigmap); /* truncate at 3 sigma-1D == sqrt(3) sigma-3D */
    g_extx = 2 * exth + 1; g_exty = g_extx; extz = g_extx; nvox = g_extx * g_exty * extz;
    
    phi = (double *) malloc(nvox * sizeof(double));
    /* write Gaussian within 3 sigma-1D to map */
    bvalue = -1.0 / (2.0*sigmap*sigmap);
    cvalue = 9.0*sigmap*sigmap;
    for(count=0;count<nvox;count++) *(phi+count) = 0.0;
    for (indz=0;indz<extz;indz++)
      for (indy=0;indy<g_exty;indy++)
	for (indx=0;indx<g_extx;indx++) {
	  dsqu = (indx-exth)*(indx-exth)+(indy-exth)*(indy-exth)+(indz-exth)*(indz-exth);
	  if (dsqu < cvalue) *(phi+g1idz(indz,indy,indx)) = kampl * exp (dsqu * bvalue);
	}
    *pg_extx = g_extx;    
    return phi;
}

void apply_kernel(
  /*input map*/ const double *phi, unsigned int nvox, unsigned int g_extx, unsigned int g_exty, unsigned int extz,
  /*kernel map*/ const double *phi2, unsigned int nvox2, unsigned int g_extx2, unsigned int g_exty2, unsigned int extz2,
  /*output map*/ double *phi3
) {
    int k,j,i,kk,jj,ii;
    int count, indv, extxy2, indz, indy, indx;
    int avgx, avgy, avgz;
    double dval;
    for (count=0;count<nvox;count++) *(phi3+count) = 0.0;
    extxy2 = g_extx2 * g_exty2;
    avgx = (g_extx2-1)/2;
    avgy = (g_exty2-1)/2;
    avgz = (extz2-1)/2;
    for (count=0;count<nvox2;count++) {
      dval = *(phi2+count);
      if (dval != 0.0) {
	indv = count;
	k = indv / extxy2;
	indv -= k * extxy2;
	j = indv / g_extx2;
	i = indv - j * g_extx2;
        /*printf("KERNEL %d %d %d %.3f\n", k-avgz,j-avgy,i-avgx,dval); */
	for (indz=0;indz<extz;indz++) {
          kk = k+indz-avgz;
          if (kk < 0 || kk >= extz) continue;
          for (indy=0;indy<g_exty;indy++) {
            jj = j+indy-avgy;
            if (jj < 0 || jj >= g_exty) continue;
            for (indx=0;indx<g_extx;indx++) { 
              ii = i+indx-avgx;
              if (ii < 0 || ii >= g_extx) continue;
	      *(phi3+g1idz(kk,jj,ii)) += *(phi+g1idz(indz,indy,indx)) * dval;
            }
          }
        }    
      }
    }
}
