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

inline void restrain_type_1(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
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
    energy += weight * 0.5 * cforce * violationsq;
    if (iab) {
      double factor = violation/dis;
      Coor force = {weight * disx * cforce*factor, 
                    weight * disy * cforce*factor,
                    weight * disz * cforce*factor};      
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
    energy += weight * 0.5 * cforce * violationsq;
    //printf("ENERGY2: %.3f\n", 0.5 * cforce * violation);
    if (iab) {
      double force = weight * cforce * violation/deff;    
      double gradfac = -1.0/6 * pow(deffsum, -7.0/6) * deff * force;
      grad_deff(gradfac, r.selection1, r.s1, r.selection2, r.s2, f);
    }
  }
}

inline void restrain_type_2(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
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
      energy += weight * 0.5 * cforce * violationsq;
    }
    else {
      double maxviol = r.par3;
      double maxviolsq = maxviol * maxviol;
      double extraviol = dis - r.par3;
      double violation2sq = maxviolsq+2*maxviol*extraviol;
      factor = maxviol/dis;
      energy += weight * 0.5 * cforce * violation2sq;
    }
    //printf("ENERGY: %.3f\n", 0.5 * cforce * violationsq);        
    if (iab) {
      Coor force = {weight * disx * cforce*factor, 
                    weight * disy * cforce*factor,
                    weight * disz * cforce*factor};      
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
      energy += weight * 0.5 * cforce * violationsq;
    }
    else {
      double maxviol = r.par3;
      double maxviolsq = maxviol * maxviol;
      double extraviol = deff - r.par3;
      double violation2sq = maxviolsq+2*maxviol*extraviol;
      factor = maxviol/deff;
      energy += weight * 0.5 * cforce * violation2sq;
    }
    
    if (iab) {
      double force = weight * cforce * factor;    
      double gradfac = -1.0/6 * pow(deffsum, -7.0/6) * deff * force;
      grad_deff(gradfac, r.selection1, r.s1, r.selection2, r.s2, f);
    }
  }
}

inline void restrain_type_3(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
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
    energy += 0.5 * weight * cforce * violationsq;
    if (iab) {
      double factor = violation/dis;
      Coor force = {weight * disx * cforce*factor, 
                    weight * disy * cforce*factor,
                    weight * disz * cforce*factor};
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

inline void restrain_type_4(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //harmonic bond restraints for preserving intramolecular secondary structure
	// poor mans elastic network force field
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection2[0]-1;
    const Coor &a2 = x[atomnr2];

    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;
    double cforce = r.par2;
    double dis = sqrt(dsq);
    double violation = dis - r.par1;
    if (dsq == 0.0) {
    	return;
    }
    double violationsq = violation*violation;
    energy += 0.5 * cforce * violationsq;
    //fprintf(stderr, "Bond restraint %i  %i %f %f \n",atomnr1, atomnr2, dis, r.par1);
    if (iab) {
    	double factor = violation/dis;
    	Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
    	//fprintf(stderr, "Bond restraint %i  %i %f %f %f %f %f\n",atomnr1, atomnr2, dis, r.par1, force[0], force[1], force[2]);
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
}

inline void restrain_type_5(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //Repulsive restraints for preserving intramolecular secondary structure
	// poor mans Lennard Jones
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
    double violation = dsq - limsq;
    double violationsq = violation*violation;
    energy += 0.5 * cforce * violationsq;
    // ToDo check force code
    if (iab) {
    	double factor = 2*violation;
    	Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
    	//fprintf(stderr, "LJ restraint %i  %i %f %f %f %f %f\n",atomnr1, atomnr2, dsq, limsq, force[0], force[1], force[2]);
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
}

inline void restrain_type_6(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //Step potential, use only for scoring and Monte Carlo!
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    int atomnr2 = r.selection2[0]-1;
    const Coor &a1 = x[atomnr1];
    const Coor &a2 = x[atomnr2];
    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;
    double limsq = r.par1 * r.par1;//upper limit of step potential
    if (limsq < dsq) {
      //printf("ENERGY: 0\n");
      return;
    }
    double limsq2 = r.par3 * r.par3;//lower limit of step potential
    if (limsq2 > dsq) {
      //printf("ENERGY: 0\n");
      return;
    }
    double cforce = r.par2;//depth of step potential
    energy += 0.5*cforce;
  }
}

inline void restrain_type_7(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //positional restraints
  Coor refe = {r.par4, r.par5, r.par6}; 
  for (int n = 0; n < r.s1; n++) {
    int atomnr1 = r.selection1[n]-1;
    const Coor &a1 = x[atomnr1];
    double dsq = 0;
    for (int i = 0; i < 3; i++) {
      char mask = (1 << i);
      if (!(r.position_type & mask)) continue;      
      double d = a1[i] - refe[i];
      dsq += d*d;
    }    
    double dminsq = r.par1 * r.par1;
    double dmaxsq = r.par2 * r.par2;
    double cforce = r.par3;
    if (dsq > dmaxsq) {      
      double dis = sqrt(dsq);
      double violation = dis - r.par2;
      double violationsq = violation*violation;
      energy += weight * 0.5 * cforce * violationsq;
      if (iab) {
        double factor = violation/dis;
        Coor &f1 = f[atomnr1];
        for (int i = 0; i < 3; i++) {
          char mask = (1 << i);
          if (!(r.position_type & mask)) continue;
          double d = a1[i] - refe[i];
          double force = weight * d * cforce * factor;
          f1[i] -= force;
        }            
      }
    }
    if (dsq < dminsq) {      
      double dis = sqrt(dsq);
      double violation = r.par1 - dis;
      double violationsq = violation*violation;
      energy += weight * 0.5 * cforce * violationsq;
      if (iab) {
        double factor = violation/dis;
        Coor &f1 = f[atomnr1];
        for (int i = 0; i < 3; i++) {
          char mask = (1 << i);
          if (!(r.position_type & mask)) continue;
          double d = a1[i] - refe[i];
          double force = weight * d * cforce * factor;
          f1[i] += force;
        }            
      }
    }    
  }
}

inline void restrain_type_8(const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //Bump potential for e.g. x-link data
  // goes smoothly to zero at r_min1 = r_ref1+r_0, r_min2= r_ref2-r_0, between r_ref1 and r_ref2 flat potential
  //can be repulsive and attracttive
  if (r.s1 == 1 && r.s2 == 1) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection2[0]-1;
    const Coor &a2 = x[atomnr2];

    double disx = a1[0]-a2[0];
    double disy = a1[1]-a2[1];
    double disz = a1[2]-a2[2];
    double dsq = disx*disx+disy*disy+disz*disz;
    double dist = sqrt(dsq);
    double limsq = (r.par1 -r.par3)*(r.par1 - r.par3);//restraint goes to zero here
    double limsq2 = (r.par2+r.par3)*(r.par2 + r.par3);//restraint goes to zero here
    double limsq3 = r.par3 * r.par3;// 
    if ((limsq2 < dsq) || (limsq > dsq) ) {
      //printf("ENERGY: 0 %f\n", dist);
      return;
    }
    limsq = r.par1 * r.par1;
    limsq2 = r.par2 * r.par2;
    double cforce = r.par4;
    
    double violation = -1.0*limsq3;
    if (dsq < limsq){
     violation += (dist - r.par1) * (dist-r.par1); 
    }
    else if (dsq > limsq2 ){
      violation += (dist - r.par2)*(dist-r.par2);
    }
    double violationsq = violation*violation;
    energy += 0.5 * cforce * violationsq;
    // ToDo check force code
    if (iab) {
    	double factor = 2*violation;
    	Coor force = {disx * cforce*factor,disy * cforce*factor,disz * cforce*factor};
	if ( dsq >= limsq && dsq <= limsq2 ) {
	  force[0] = 0.0;
	  force[1] = 0.0;
	  force[2] = 0.0;
	}
	  //fprintf(stderr, "LJ restraint %i  %i %f %f %f %f %f\n",atomnr1, atomnr2, dsq, limsq, force[0], force[1], force[2]);
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
}
inline void restrain_type_9(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //ANGLE-type restaint between 3 atoms
  if (r.s1 == 3 ) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection1[1]-1;
    const Coor &a2 = x[atomnr2];
    int atomnr3 = r.selection1[2]-1;
    const Coor &a3 = x[atomnr3];
    double x12= a1[0]-a2[0];
    double y12= a1[1]-a2[1];
    double z12= a1[2]-a2[2];
    double x32= a3[0]-a2[0];
    double y32= a3[1]-a2[1];
    double z32= a3[2]-a2[2];
    double  dis2 = 1.0/(x12*x12+y12*y12+z12*z12);
    double  dis = sqrt(dis2);
    double  des2 = 1.0/(x32*x32+y32*y32+z32*z32);
    double  des = sqrt(des2);
    double  skalar = x12*x32+y12*y32+z12*z32;
    double  winkel = acos(skalar*des*dis);
    double  diff = winkel-r.par1;
    double  cforce = r.par2;
    double  fdiff = cforce*diff;
    energy += weight *  fdiff * diff;
//   printf("W %.3f %3f %3f %3f \n", winkel,cforce,r.par1,r.par2);
    if (iab) {
      diff=dis*des;
      dis=skalar*diff;
      double factor = -2.0*fdiff/sqrt(1.0 - dis*dis);
      diff=factor*diff;
      double diffdis=factor*dis*dis2;
      Coor force1 = {weight * (x32*diff-diffdis*x12),
                    weight * (y32*diff-diffdis*y12),
                    weight * (z32*diff-diffdis*z12)};
      Coor &f1 = f[atomnr1];
      f1[0] -= force1[0];
      f1[1] -= force1[1];
      f1[2] -= force1[2];
      Coor &f2 = f[atomnr2];
      f2[0] += force1[0];
      f2[1] += force1[1];
      f2[2] += force1[2];
      double diffdes=factor*dis*des2;
      Coor force2 = {weight * (x12*diff-diffdes*x32),
                    weight * (y12*diff-diffdes*y32),
                    weight * (z12*diff-diffdes*z32)};
      Coor &f3 = f[atomnr3];
      f3[0] -= force2[0];
      f3[1] -= force2[1];
      f3[2] -= force2[2];
      f2[0] += force2[0];
      f2[1] += force2[1];
      f2[2] += force2[2];
    }
  }
}
inline void restrain_type_10(double weight, const Restraint &r, int iab, const Coor *x, Coor *f, double &energy) {
  //Dihedral-type restaint between 4 atoms
  if (r.s1 == 4 ) {
    int atomnr1 = r.selection1[0]-1;
    const Coor &a1 = x[atomnr1];
    int atomnr2 = r.selection1[1]-1;
    const Coor &a2 = x[atomnr2];
    int atomnr3 = r.selection1[2]-1;
    const Coor &a3 = x[atomnr3];
    int atomnr4 = r.selection1[3]-1;
    const Coor &a4 = x[atomnr4];
    double x21= a2[0]-a1[0];
    double y21= a2[1]-a1[1];
    double z21= a2[2]-a1[2];
    double x32= a2[0]-a3[0];
    double y32= a2[1]-a3[1];
    double z32= a2[2]-a3[2];
    double x31= a3[0]-a1[0];
    double y31= a3[1]-a1[1];
    double z31= a3[2]-a1[2];
    double x43= a4[0]-a3[0];
    double y43= a4[1]-a3[1];
    double z43= a4[2]-a3[2]; 
    double x42= a4[0]-a2[0];
    double y42= a4[1]-a2[1];
    double z42= a4[2]-a2[2];
    double r13x=z21*y32-y21*z32;
    double r13y=x21*z32-z21*x32;
    double r13z=y21*x32-x21*y32;
    double r24x=z32*y43-y32*z43;
    double r24y=x32*z43-z32*x43;
    double r24z=y32*x43-x32*y43;
    double  r13n = 1.0/sqrt(r13x*r13x+r13y*r13y+r13z*r13z);
    double  r24n = 1.0/sqrt(r24x*r24x+r24y*r24y+r24z*r24z);
    r13x=r13x*r13n;
    r13y=r13y*r13n;
    r13z=r13z*r13n;
    r24x=r24x*r24n;
    r24y=r24y*r24n;
    r24z=r24z*r24n;
    double cijkl=r13x*r24x+r13y*r24y+r13z*r24z;
    if (cijkl >= 1.0) cijkl=0.999999999;
    if (cijkl < -1.0) cijkl=-0.999999999;
    double sijkl=x32*(r13z*r24y-r13y*r24z)+y32*(r13x*r24z-r13z*r24x)+z32*(r13y*r24x-r13x*r24y);
    if (sijkl < 0 ) sijkl=-1.0; else sijkl=1.0;
    double  winkel = sijkl*acos(cijkl);
    if (winkel >= 0 && winkel < 0.000001) winkel = 0.000001;
    if (winkel < 0 && winkel > -0.000001) winkel = -0.000001;
    double  delwink = winkel-r.par1;      
    if(delwink > 3.141592654) delwink=delwink-6.283185307;
    if(delwink < -3.141592654) delwink=delwink+6.283185307;
    double factor=r.par2*delwink;
    energy += weight * factor*delwink;
//    printf("W %.3f %3f %3f %3f \n", winkel,r.par1,r.par2,delwink);
    if (iab) {
      factor=2.0 * weight * factor/sin(winkel);
      double v1x=r13x-cijkl*r24x;
      double v1y=r13y-cijkl*r24y;
      double v1z=r13z-cijkl*r24z;
      double v2x=r24x-cijkl*r13x;
      double v2y=r24y-cijkl*r13y;
      double v2z=r24z-cijkl*r13z;
      double v3x=z32*v2y-y32*v2z;
      double v3y=x32*v2z-z32*v2x;
      double v3z=y32*v2x-x32*v2y;
      Coor fi = {factor * r13n * v3x,
                 factor * r13n * v3y,
                 factor * r13n * v3z};
      v3x=z31*v2y-y31*v2z;
      v3y=x31*v2z-z31*v2x;
      v3z=y31*v2x-x31*v2y;
      double v4x=z43*v1y-y43*v1z;
      double v4y=x43*v1z-z43*v1x;
      double v4z=y43*v1x-x43*v1y;
      Coor fj = {factor * ( r13n * v3x - r24n * v4x),
                 factor * ( r13n * v3y - r24n * v4y),
                 factor * ( r13n * v3z - r24n * v4z)};
      v3x=z42*v1y-y42*v1z;
      v3y=x42*v1z-z42*v1x;
      v3z=y42*v1x-x42*v1y;
      v4x=z21*v2y-y21*v2z;
      v4y=x21*v2z-z21*v2x;
      v4z=y21*v2x-x21*v2y;
      Coor fk = {factor * ( r24n * v3x - r13n * v4x),
                 factor * ( r24n * v3y - r13n * v4y),
                 factor * ( r24n * v3z - r13n * v4z)};
      v3x=z32*v1y-y32*v1z;
      v3y=x32*v1z-z32*v1x;
      v3z=y32*v1x-x32*v1y;
      Coor fl = {factor *  r24n * v3x,
                factor *  r24n * v3y,
                factor *  r24n * v3z };
      Coor &f1 = f[atomnr1];
      f1[0] -= fi[0];
      f1[1] -= fi[1];
      f1[2] -= fi[2];
      Coor &f2 = f[atomnr2];
      f2[0] -= fj[0];
      f2[1] -= fj[1];
      f2[2] -= fj[2];
      Coor &f3 = f[atomnr3];
      f3[0] -= fk[0];
      f3[1] -= fk[1];
      f3[2] -= fk[2];
      Coor &f4 = f[atomnr4];
      f4[0] -= fl[0];
      f4[1] -= fl[1];
      f4[2] -= fl[2];
    }
  }
}

extern "C" void restrain_(const int &ministatehandle, const int &cartstatehandle, const int &seed, const int &iab, double &energy) {
  MiniState &ms = ministate_get(ministatehandle);
  CartState &cs = cartstate_get(cartstatehandle);
  Coor *x = (Coor *) &cs.x[0];
  Coor *f = (Coor *) &cs.f[0];
  double weight = ms.restweight;
  srand(seed);
  for (int n = 0; n < ms.nr_restraints; n++) {
    Restraint &r = ms.restraints[n];
    if (r.maxindex > cs.nall) {
      fprintf(stderr, "Restraint %d selects atoms up to index %d, but there are only %d atoms\n", n+1, r.maxindex, cs.nall);
      exit(1);
    }
    if (r.type == 1) restrain_type_1(weight, r,iab,x,f,energy);
    if (r.type == 2) restrain_type_2(weight,r,iab,x,f,energy);    
    if (r.type == 3) restrain_type_3(weight,r,iab,x,f,energy);    
    if (r.type == 4) restrain_type_4(r,iab,x,f,energy);
    if (r.type == 5) restrain_type_5(r,iab,x,f,energy);
    if (r.type == 6) restrain_type_6(r,iab,x,f,energy);
    if (r.type == 7) restrain_type_7(weight, r,iab,x,f,energy);
    if (r.type == 8) restrain_type_8(r,iab,x,f,energy);
    if (r.type == 9) restrain_type_9(weight,r,iab,x,f,energy);
    if (r.type == 10) restrain_type_10(weight,r,iab,x,f,energy);  
  }
}
