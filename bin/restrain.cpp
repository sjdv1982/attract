#include "state.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

extern MiniState &ministate_get(int handle);
extern CartState &cartstate_get(int handle);

struct DeffStruct {
  double f;
  Coor d;
};

static DeffStruct deffs[MAXSELECTION*MAXSELECTION];

inline double calc_deffsum(int iab, const int *selection1, int s1, const int *selection2, int s2, const Coor *x, double lim) {

  bool has_lim = 0;
  double deffsum_lim = 0;
  if (lim > 0) {
    has_lim = 1;
    double limsq = lim*lim;
    deffsum_lim = 1.0/(limsq*limsq*limsq);
  }
  double deffsum = 0;
  
  DeffStruct *curr_deff = deffs;
  for (int n1 = 0; n1 < s1; n1++) {
    int atomnr1 = selection1[n1]-1;
    const Coor &a1 = x[atomnr1];
    for (int n2 = 0; n2 < s2; n2++) {
      int atomnr2 = selection2[n2]-1;
      const Coor &a2 = x[atomnr2];
      double disx = a1[0]-a2[0];
      double disy = a1[1]-a2[1];
      double disz = a1[2]-a2[2];
      double dsq = disx*disx+disy*disy+disz*disz;      
      double deffsum_contrib = 1.0/(dsq * dsq * dsq);
      deffsum += deffsum_contrib;
      if (has_lim) if (deffsum >= deffsum_lim) return -1;
      if (iab) {
        curr_deff->f = dsq;
        curr_deff->d[0] = disx; 
        curr_deff->d[1] = disy; 
        curr_deff->d[2] = disz;       
        curr_deff++;
      }
    }
  } 
  
  if (iab) {
    DeffStruct *curr_deff = deffs;
    for (int n1 = 0; n1 < s1; n1++) {
      for (int n2 = 0; n2 < s2; n2++) {
	double dsq = curr_deff->f;
	double ddsq = 1.0/dsq;
	//derivative for dsq: -2 * d
	//derivative for dsq^-3: -3 * dsq^-4 * -2 * d = 6 * dsq^-4 * d
	//divide by d (will be later multiplied)
	curr_deff->f = 6*ddsq*ddsq*ddsq*ddsq;
	
	curr_deff++;
      }
    }  
  }
  return deffsum;
} 

inline void grad_deff(double gradfac, const int *selection1, int s1, const int *selection2, int s2, Coor *f) {

  DeffStruct *curr_deff = deffs;
  for (int n1 = 0; n1 < s1; n1++) {
    int atomnr1 = selection1[n1]-1;
    Coor &a1 = f[atomnr1];
    for (int n2 = 0; n2 < s2; n2++) {
      int atomnr2 = selection2[n2]-1;
      Coor &a2 = f[atomnr2];
      double disfac = curr_deff->f * gradfac;
      double gradx = curr_deff->d[0] * disfac; 
      double grady = curr_deff->d[1] * disfac; 
      double gradz = curr_deff->d[2] * disfac;             
      a1[0] += gradx; a1[1] += grady; a1[2] += gradz;
      a2[0] -= gradx; a2[1] -= grady; a2[2] -= gradz;
      curr_deff++;
    }
  }  

}

inline void restrain_type_1(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //maximum-distance harmonic restraints
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection2[0]-1;
    const Coor &a2 = x[atomnr2];

    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;      
    double limsq = r.par1 * r.par1;
    if (limsq > dsq) {
      //printf("ENERGY: 0\n");
      return;
    }
    double cforce = r.par2;
    double dis = sqrt(dsq);
    double violation = dis - r.par1;
    double violationsq = violation*violation;
    energy += 0.5 * cforce * violationsq;
    if (iab) {
      double factor = violation/dis;
      Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
      Coor &f1 = f[atomnr1];
      f1[0] -= force[0];
      f1[1] -= force[1];
      f1[2] -= force[2];      
      Coor &f2 = f[atomnr2];
      f2[0] += force[0];
      f2[1] += force[1];
      f2[2] += force[2];
    }
  }
  else {
    double deffsum = calc_deffsum(iab, r.selection1, r.s1, r.selection2, r.s2, x, r.par1);
    if (deffsum == -1) {
      //printf("ENERGY2: 0\n");
      return;
    }
    double deff = pow(deffsum, -1.0/6);
    double violation = deff - r.par1;
    double violationsq = violation*violation;
    double cforce = r.par2;
    energy += 0.5 * cforce * violationsq;
    //printf("ENERGY2: %.3f\n", 0.5 * cforce * violation);
    if (iab) {
      double force = cforce * violation/deff;    
      double gradfac = -1.0/6 * pow(deffsum, -7.0/6) * deff * force;
      grad_deff(gradfac, r.selection1, r.s1, r.selection2, r.s2, f);
    }
  }
}

inline void restrain_type_2(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //HADDOCK-type AIR restraints
  if (double(rand())/RAND_MAX <= r.par4) return;
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection2[0]-1;
    const Coor &a2 = x[atomnr2];

    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;      
    double limsq = r.par1 * r.par1;
    if (limsq > dsq) {
      //printf("ENERGY: 0\n");
      return;
    }
    double cforce = r.par2;
    double dis = sqrt(dsq);
    double violation = dis - r.par1;
    double violationsq = violation*violation;
    double factor;
    if (violation < r.par3) {
      factor = violation/dis;
      energy += 0.5 * cforce * violationsq;
    }
    else {
      double maxviol = r.par3;
      double maxviolsq = maxviol * maxviol;
      double extraviol = dis - r.par3;
      double violation2sq = maxviolsq+2*maxviol*extraviol;
      factor = maxviol/dis;
      energy += 0.5 * cforce * violation2sq;
    }
    //printf("ENERGY: %.3f\n", 0.5 * cforce * violationsq);        
    if (iab) {
      Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
      Coor &f1 = f[atomnr1];
      f1[0] -= force[0];
      f1[1] -= force[1];
      f1[2] -= force[2];      
      Coor &f2 = f[atomnr2];
      f2[0] += force[0];
      f2[1] += force[1];
      f2[2] += force[2];
    }
  }
  else {
    double deffsum = calc_deffsum(iab, r.selection1, r.s1, r.selection2, r.s2, x, r.par1);
    if (deffsum == -1) {
      //printf("ENERGY2: 0\n");
      return;
    }
    double deff = pow(deffsum, -1.0/6);
    
    double cforce = r.par2;
    double violation = deff - r.par1;
    double violationsq = violation*violation;
    double factor;
    if (violation < r.par3) {
      factor = violation/deff;
      energy += 0.5 * cforce * violationsq;
    }
    else {
      double maxviol = r.par3;
      double maxviolsq = maxviol * maxviol;
      double extraviol = deff - r.par3;
      double violation2sq = maxviolsq+2*maxviol*extraviol;
      factor = maxviol/deff;
      energy += 0.5 * cforce * violation2sq;
    }
    
    if (iab) {
      double force = cforce * factor;    
      double gradfac = -1.0/6 * pow(deffsum, -7.0/6) * deff * force;
      grad_deff(gradfac, r.selection1, r.s1, r.selection2, r.s2, f);
    }
  }
}

inline void restrain_type_3(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //minimum-distance harmonic restraints
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection2[0]-1;
    const Coor &a2 = x[atomnr2];

    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;      
    double limsq = r.par1 * r.par1;
    if (limsq < dsq) {
      //printf("ENERGY: 0\n");
      return;
    }
    double cforce = r.par2;
    double dis = sqrt(dsq);
    double violation = r.par1 - dis;
    double violationsq = violation*violation;
    energy += 0.5 * cforce * violationsq;
    if (iab) {
      double factor = violation/dis;
      Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
      Coor &f1 = f[atomnr1];
      f1[0] += force[0];
      f1[1] += force[1];
      f1[2] += force[2];      
      Coor &f2 = f[atomnr2];
      f2[0] -= force[0];
      f2[1] -= force[1];
      f2[2] -= force[2];
    }
  }
}


extern "C" void restrain_(const int &ministatehandle, const int &cartstatehandle, const int &seed, const int &iab, double &energy) {
  MiniState &ms = ministate_get(ministatehandle);
  CartState &cs = cartstate_get(cartstatehandle);
  Coor *x = (Coor *) &cs.x[0];
  Coor *f = (Coor *) &cs.f[0];
  srand(seed);
  for (int n = 0; n < ms.nr_restraints; n++) {
    Restraint &r = ms.restraints[n];
    if (r.type == 1) restrain_type_1(r,iab,x,f,energy);
    if (r.type == 2) restrain_type_2(r,iab,x,f,energy);    
    if (r.type == 3) restrain_type_3(r,iab,x,f,energy);    
  }
}
