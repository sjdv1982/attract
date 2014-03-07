#include "axsym.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

const double pi = 4.0f * atan(1.0f);

extern void axsym_fold_grads(
 MiniState &m, CartState &c,
 double *grads0, double *grads, double *morphgrads
) 
{
  if (c.nr_symtrans == 0) {
    memcpy(grads, grads0, MAXDOF * sizeof(double));
    return;
  }
  int jb0 = 3*m.iori*(c.nlig0-m.fixre)+3*m.itra*(c.nlig0-m.fixre);
  int jb = 3*m.iori*(c.nlig-m.fixre)+3*m.itra*(c.nlig-m.fixre);
  int jl0=3*m.iori*(c.nlig0-m.fixre);
  int jl=3*m.iori*(c.nlig-m.fixre);
  int ju = jb;
  int ju0 = jb0;
  for (int idl1 =0; idl1 < c.nlig0; idl1++){
	  ju0 += c.nhm[idl1];
  }
  for (int idl1 =0; idl1 < c.nlig; idl1++){
  	  ju += c.nhm[idl1];
    }
  //printf("NLIG0 JB0 JL0 IORI ITRA FIXRE %d %d %d %d %d %d\n", c.nlig0, jb0, jl0, m.iori, m.itra, m.fixre);
  //printf("NLIG JB JL IORI ITRA FIXRE %d %d %d %d %d %d\n", c.nlig, jb, jl, m.iori, m.itra, m.fixre);
  for (int idl1 = 0; idl1 < c.nlig0; idl1++) {
    for (int cop = 0; cop < c.nr_symcopies[idl1]; cop++) {
      int idl2 = c.symcopies[idl1][cop];
      //printf("IDL1 IDL2 %d %d\n", idl1, idl2);
      int pos1, pos2;
      if (m.iori && idl1 >= m.fixre) {
        pos1 = 3 * (idl1-m.fixre);
        pos2 = 3 * (idl2-m.fixre);
        //printf("ORI POS1 POS2 %d %d %.3f\n", pos1, pos2, grads0[pos2]);
        grads[pos1] += grads0[pos2];
        grads[pos1+1] += grads0[pos2+1];
        grads[pos1+2] += grads0[pos2+2];
      }

      if (m.itra && idl1 >= m.fixre) {
        pos1 = jl0 + 3 * (idl1-m.fixre);
        pos2 = jl + 3 * (idl2-m.fixre);
        //printf("TRANS POS1 POS2 %d %d %.3f\n", pos1, pos2, grads0[pos2]);
        grads[pos1] += grads0[pos2];
        grads[pos1+1] += grads0[pos2+1];
        grads[pos1+2] += grads0[pos2+2];
      }
      if (m.ieig) {
        int pos1 = jb0;
        for (int nn = 0; nn < idl1; nn++) {
          pos1 += c.nhm[nn];
        }
        int pos2 = jb;
        for (int nn = 0; nn < idl2; nn++) {
          pos2 += c.nhm[nn];
        }
        for (int nn = 0; nn < c.nhm[idl1]; nn++) {
          grads[pos1+nn] = grads0[pos2+nn];
        }
      }
      if (m.iindex) {
    	  int pos1 = ju0;
    	  for (int nn = 0; nn < idl1; nn++) {
    	      pos1 += c.nihm[nn];
    	  }
    	  int pos2 = ju;
    	  for (int nn = 0; nn < idl2; nn++) {
    	      pos2 += c.nihm[nn];
    	   }
    	   for (int nn = 0; nn < c.nihm[idl1]; nn++) {
    	      grads[pos1+nn] = grads0[pos2+nn];
    	   }
      }

      morphgrads[idl1] = morphgrads[idl2];
    }  
  }
}
extern void prepare_axsym_cartstate(CartState &c) { 
  int nlig_dof = prepare_axsym_dof(
   c.nr_axsyms,
   c.axsyms,
   c.nlig0,
   c.nhm,
   c.nihm,
   c.nrens,
   c.morphing,
   c.has_locrests,

   c.nr_symcopies,
   c.symcopies,
   c.nr_symtrans,
   c.symtrans,
   
   c.forcefactor,
   c.forcerotation
  );
  int nlig = c.nlig0;
  for (int n = 0; n < c.nr_symtrans; n++) {
    SymTrans &sym = c.symtrans[n];
    int l = sym.ligand;
    int target = sym.targetligand;
    if (l == target) continue;
    if (target != nlig) {
      fprintf(stderr, "AssertionError axsym.cpp, target: %d, nlig: %d\n", target, nlig);
      exit(1);
    }
    
    int natoml = c.natom[l];
    int n3atoml = c.n3atom[l];
    int nresl = c.nres[l];
    int pos = 0;
    if (l > 0) pos = c.ieins[l-1];
    int pos3 = 0;
    if (l > 0) pos3 = c.ieins3[l-1];
    int posres = 0;
    for (int nn = 0; nn < l; nn++) posres += c.nres[nn];

    int cnall = c.nall;
    int cnall3 = c.nall3;
    int nresall = 0;    
    
    memcpy(&c.iei[cnall], &c.iei[pos], natoml*sizeof(int)); 
    memcpy(&c.x[cnall3], &c.x[pos3], n3atoml*sizeof(double)); 
    c.pivot[0][nlig] = c.pivot[0][l];
    c.pivot[1][nlig] = c.pivot[1][l];
    c.pivot[2][nlig] = c.pivot[2][l];
    memcpy(&c.xb[cnall3], &c.xb[pos3], n3atoml*sizeof(double)); 
    memcpy(&c.xori0[cnall3], &c.xori0[pos3], n3atoml*sizeof(double)); 
    memcpy(&c.xori[cnall3], &c.xori[pos3], n3atoml*sizeof(double)); 
    memcpy(&c.iaci[cnall], &c.iaci[pos], natoml*sizeof(int)); 
    memcpy(&c.iaci_old[cnall], &c.iaci_old[pos], natoml*sizeof(int)); 
    memcpy(&c.xlai[cnall], &c.xlai[pos], natoml*sizeof(double)); 
    memcpy(&c.icop[cnall], &c.icop[pos], natoml*sizeof(int)); 
    memcpy(&c.we[cnall], &c.we[pos], natoml*sizeof(double)); 
    memcpy(&c.chai[cnall], &c.chai[pos], natoml*sizeof(double)); 

    for (int i = 0; i < c.nhm[l]; i++) {
      c.val[i][nlig] = c.val[i][l];
    }
    for (int j = 0; j < n3atoml; j++) {      
      for (int i = 0; i < c.nhm[l]; i++) {
        c.eig[j][i][nlig] = c.eig[j][i][l];
      }
    }
    for (int i = 0; i < c.nihm[l]; i++) {
    	for (int j=0; j < MAXLENINDEXMODE; j++){
    		c.index_eig[j][i][nlig] = c.index_eig[j][i][l];
    		c.index_val[j][i][nlig] = c.index_val[j][i][l];
    	}
    }
    memcpy(&c.ncop[nresall], &c.ncop[posres], nresl*21*11*sizeof(int)); 
    memcpy(&c.nmaxco[nresall], &c.nmaxco[posres], nresl*sizeof(int)); 
    memcpy(&c.natco[nresall], &c.natco[posres], nresl*sizeof(int)); 
    c.grids[nlig] = c.grids[l];
    memcpy(c.ensd[nlig], c.ensd[l], c.nrens[l] * sizeof(double *));
    memcpy(c.morphd[nlig], c.morphd[l], c.nrens[l] * sizeof(double *));
    c.ensw[nlig] = c.ensw[l];
    cnall += natoml;
    cnall3 += n3atoml;
    nresall += nresl;
    c.natom[nlig] = natoml;
    c.n3atom[nlig] = n3atoml;
    c.ieins[nlig] = cnall;
    c.ieins3[nlig] = cnall3;                
    nlig++;

    c.nall = cnall;
    c.nall3 = cnall3;
  }
  c.nlig = nlig;
  if (nlig != nlig_dof) {
    fprintf(stderr, "AssertionError axsym.cpp, nlig != nlig_dof: %d vs %d\n", nlig, nlig_dof);
    exit(1);
  };
}

int prepare_axsym_dof(
 int nr_axsyms,
 AxSymmetry *axsyms,
 int nlig0,

 int *nhm,
 int *nihm,
 int *nrens,
 int *morphing,
 int *has_locrests,
 
 int *nr_symcopies,
 int (&symcopies)[MAXLIG][24*MAXLIG],
 int &nr_symtrans,
 SymTrans *symtrans,
 
 double *forcefactor,
 double (&forcerotation)[MAXLIG][9]
)
{
  nr_symtrans = 0;
  int nlig = nlig0;
  for (int n = 0; n < nr_axsyms; n++) {
    AxSymmetry &sym = axsyms[n];
    int nrsym = sym.symtype;
    int l = sym.ligand - 1;
    if (nrsym < 0) {
      fprintf(stderr, "Symmetry type must be positive!\n");
      exit(1);
    }
    int cnr_symcopies = nr_symcopies[l];
    for (int cop = 0; cop < cnr_symcopies; cop++) {
      int ll = symcopies[l][cop];
      if (forcefactor[ll] == 0) forcefactor[ll] = 1;
      double ff = forcefactor[ll] * nrsym;
      forcefactor[ll] = ff;
      int symcount = nrsym;
      if (nrsym == 0) symcount = 2; //ncsym => just run the loop once
      for (int nn = 1; nn < symcount; nn++) {

        //Copy DOF descriptors from the original
        nhm[nlig] = nhm[l];
        nihm[nlig] = nihm[l];
        nrens[nlig] = nrens[l]; 
        morphing[nlig] = morphing[l],      
        has_locrests[nlig] = has_locrests[l];            
        
        //Create a new symmetry transformation
        SymTrans &csymtrans = symtrans[nr_symtrans];
        csymtrans.ligand = ll;
        csymtrans.targetligand = nlig;
        memcpy(csymtrans.origin, sym.origin, 3 * sizeof(double));
        
        double angle;
        if (nrsym == 0) { //ncsym
          angle = sym.angle / 180 * pi;
        }
        else {
          angle = 2.0 * pi / nrsym * nn;        
        }

        //Compute a symmetry matrix from the symmetry axis and the current angle
        double *rotmatsym = csymtrans.rotmatsym;
        double c = cos(angle);
        double s = sin(angle);
        double t = 1 - c;
        double x = sym.axis[0];
        double y = sym.axis[1];
        double z = sym.axis[2];
        rotmatsym[0] = t*x*x + c;
        rotmatsym[1] = t*x*y + z*s;
        rotmatsym[2] = t*x*z - y*s;
        rotmatsym[3] = t*x*y - z*s;
        rotmatsym[4] = t*y*y + c;
        rotmatsym[5] = t*y*z + x*s;
        rotmatsym[6] = t*x*z + y*s;
        rotmatsym[7] = t*y*z - x*s;
        rotmatsym[8] = t*z*z + c; 
	//printf("%.3f %.3f %.3f\n", rotmatsym[0], rotmatsym[1], rotmatsym[2]);
	//printf("%.3f %.3f %.3f\n", rotmatsym[3], rotmatsym[4], rotmatsym[5]);
	//printf("%.3f %.3f %.3f\n\n", rotmatsym[6], rotmatsym[7], rotmatsym[8]);

        nr_symtrans++;
        symcopies[l][nr_symcopies[l]] = nlig;        
        nr_symcopies[l]++;

        forcefactor[nlig] = ff;
        matmult_(forcerotation[ll], rotmatsym, forcerotation[nlig]);
        
        nlig++;
      }
    }  
  } 
  return nlig;
}

void apply_axsym(
 int nr_symtrans,
 SymTrans *symtrans,

 double *morph, 
 int *ens, 
 double *phi,
 double *ssi,
 double *rot,
 double *xa,
 double *ya,
 double *za,
 modes2 &dlig, 
 const int *has_locrests,
 coors2 &locrests
) 
{
  if (nr_symtrans == 0) return;
  
  for (int n = 0; n < nr_symtrans; n++) {
    SymTrans &csymtrans = symtrans[n];
    int l = csymtrans.ligand;
    int target = csymtrans.targetligand;
    double *rotmatsym = csymtrans.rotmatsym;
    double *origin = csymtrans.origin;
    
    //Compute ligand rotation matrix
    double rotmatl[9];
    euler2rotmat_(phi[l],ssi[l],rot[l],rotmatl);    
            
    //Rotate the center using the symmetry matrix
    double lc[3];
    lc[0] = xa[l] - origin[0];
    lc[1] = ya[l] - origin[1];
    lc[2] = za[l] - origin[2];
    double lc2[3];
    vecmatmult_(lc, rotmatsym,lc2);
    xa[target] = lc2[0] + origin[0];
    ya[target] = lc2[1] + origin[1];
    za[target] = lc2[2] + origin[2];

    //Multiply the ligand's original rotation matrix with the symmetry matrix
    double rotmatd[9];
    matmult_(rotmatl, rotmatsym,rotmatd);

    
    //printf("%.3f %.3f %.3f\n", rotmatd[0], rotmatd[1], rotmatd[2]);
    //printf("%.3f %.3f %.3f\n", rotmatd[3], rotmatd[4], rotmatd[5]);
    //printf("%.3f %.3f %.3f\n\n", rotmatd[6], rotmatd[7], rotmatd[8]);
    
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

    //Copy other DOFs from the original
    morph[target] = morph[l];
    memcpy(dlig[target], dlig[l], (MAXMODE+MAXINDEXMODE)*sizeof(double));
    if (has_locrests[l]) {
      memcpy(locrests[target], locrests[l], 3*sizeof(double)); 
    }  
  }
  return;  
}
