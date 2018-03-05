#include "math.h"
#define pi           3.14159265358979323846  /* pi */
#include "axsym.h"

typedef double Rotmat[9];
typedef double Vec[3];

inline void matmult(const Rotmat &rotmat1, const Rotmat &rotmat2, Rotmat &rotmat) {
   const Rotmat &m1 = rotmat1; const Rotmat &m2 = rotmat2; Rotmat &m = rotmat;
   m[0] = m1[0]*m2[0]+m1[3]*m2[1]+m1[6]*m2[2];
   m[1] = m1[1]*m2[0]+m1[4]*m2[1]+m1[7]*m2[2];
   m[2] = m1[2]*m2[0]+m1[5]*m2[1]+m1[8]*m2[2];
   m[3] = m1[0]*m2[3]+m1[3]*m2[4]+m1[6]*m2[5];
   m[4] = m1[1]*m2[3]+m1[4]*m2[4]+m1[7]*m2[5];
   m[5] = m1[2]*m2[3]+m1[5]*m2[4]+m1[8]*m2[5];
   m[6] = m1[0]*m2[6]+m1[3]*m2[7]+m1[6]*m2[8];
   m[7] = m1[1]*m2[6]+m1[4]*m2[7]+m1[7]*m2[8];
   m[8] = m1[2]*m2[6]+m1[5]*m2[7]+m1[8]*m2[8];
}

inline void euler2rotmat(double phi,double ssi,double rot, Rotmat &rotmat) {
  double cs = cos(ssi);
  double cp = cos(phi);
  double ss = sin(ssi);
  double sp = sin(phi);
  double cscp = cs*cp;
  double cssp = cs*sp;
  double sscp = ss*cp;
  double sssp = ss*sp;
  double crot = cos(rot);
  double srot = sin(rot);
  rotmat[0] = crot * cscp + srot * sp;  
  rotmat[1] = srot * cscp - crot * sp;
  rotmat[2] = sscp;
  
  rotmat[3] = crot * cssp - srot * cp;
  rotmat[4] = srot * cssp + crot * cp;
  rotmat[5] = sssp;
  
  rotmat[6] = -crot * ss;
  rotmat[7] = -srot * ss;
  rotmat[8] = cs;  
}

inline void vecmatmult(Vec &v0, Rotmat &m, Vec &v) {
  v[0] = v0[0]*m[0]+v0[1]*m[1]+v0[2]*m[2];
  v[1] = v0[0]*m[3]+v0[1]*m[4]+v0[2]*m[5];
  v[2] = v0[0]*m[6]+v0[1]*m[7]+v0[2]*m[8];
}

extern "C" void apply_axsym(
 int nr_symtrans,
 SymTrans *symtrans,
 int nstruc,
 int nlig,
 double *all_phi,
 double *all_ssi,
 double *all_rot,
 double *all_xa,
 double *all_ya,
 double *all_za
) 
{
  for (int n = 0; n < nr_symtrans; n++) {
    SymTrans &csymtrans = symtrans[n];
    int l = csymtrans.ligand;
    int target = csymtrans.targetligand;
    Rotmat &rotmatsym = csymtrans.rotmatsym;
    Vec &origin = csymtrans.origin;
    
    for (int struc = 0; struc < nstruc; struc++) {
      int offset = nlig * struc;
      double *phi = all_phi + offset;
      double *ssi = all_ssi + offset;
      double *rot = all_rot + offset;
      double *xa = all_xa + offset;
      double *ya = all_ya + offset;
      double *za = all_za + offset;
      
      //Compute ligand rotation matrix
      Rotmat rotmatl;
      euler2rotmat(phi[l],ssi[l],rot[l],rotmatl);    
              
      //Rotate the center using the symmetry matrix
      Vec lc;
      lc[0] = xa[l] - origin[0];
      lc[1] = ya[l] - origin[1];
      lc[2] = za[l] - origin[2];
      Vec lc2;
      vecmatmult(lc, rotmatsym,lc2);
      xa[target] = lc2[0] + origin[0];
      ya[target] = lc2[1] + origin[1];
      za[target] = lc2[2] + origin[2];

      //Multiply the ligand's original rotation matrix with the symmetry matrix
      double rotmatd[9];
      matmult(rotmatl, rotmatsym,rotmatd);
      
      //Distill the euler angles from the new matrix
      phi[target] = atan2(rotmatd[5],rotmatd[2]);
      ssi[target] = acos(rotmatd[8]);
      rot[target] = atan2(-rotmatd[7],-rotmatd[6]);       
      if (fabs(rotmatd[6]) < 0.0001) rot[target] = pi;
      if (fabs(rotmatd[8]) >= 0.9999) { //gimbal lock
        phi[target] = 0;
        if (fabs(rotmatd[0]) >= 0.9999) {
          if (rotmatd[0] * rotmatd[8] < 0) {
            rot[target] = pi;
          }
          else {
            rot[target] = 0;      
          }        
          if (rotmatd[8] < 0) {
            ssi[target] = pi;     
          }
          else {
            ssi[target] = 0;      
          }        
        }
        else {
          if (rotmatd[8] < 0) {
            ssi[target] = pi;     
            rot[target] = -acos(-rotmatd[0]);
          }
          else {
            ssi[target] = 0;      
            rot[target] = acos(rotmatd[0]);
          }
        }
        if (rotmatd[1] < 0) rot[target] *= -1;        
      }
    }
  }
  return;  
}
