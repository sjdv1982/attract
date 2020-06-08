#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "state.h"
#include "axsym.h"

extern "C" void read_parameters_(char *paramfile, double *rbc, 
 double *rc,double *ac,double *emin,double *rmin2,int *ipon,int *haspar, int &potshape, float &swi_on, float &swi_off,
 int paramfile_len);

CartState *cartstates[100];
int cartstatesize = 0;

extern bool exists(const char *f);

extern "C" void cartstate_get_parameters_(const int &handle,double *&rbc, 
  double *&rc,double *&ac,double *&emin,double *&rmin2,int *&ipon, int &potshape, int &cdie, double &epsilon, float &swi_on, float &swi_off, int *&use_softcore, double *&softcore);

int cartstate_new(int argc, char *argv[], bool single=0) {
  CartState *s0 = new CartState;
  memset(s0,0,sizeof(CartState));
  CartState &s = *s0;
  cartstates[cartstatesize] = s0;
  cartstatesize++;
  int cartstatehandle = cartstatesize-1+9990; 
  memset(s.haspar,0,sizeof(iParameters));
  s.transtable = NULL;
  s.epsilon = 15; 
  s.cdie = 0;
  s.morph_fconstant = 1;
  s.use_softcore = 0; 
  s.softcore=0; 
  int i;
  int dmmy; float dmmy2, dmmy3; double dmmy4;
  if (argv[0] != NULL) {
    double *rbc; double *rc; double *ac; int *use_softcore; double *softcore;
    double *emin; double *rmin2; int *ipon;
    cartstate_get_parameters_(cartstatehandle,
      rbc,rc,ac,emin,rmin2,ipon,dmmy, dmmy, dmmy4, dmmy2, dmmy3, use_softcore, softcore);
read_parameters_(argv[0],rbc,rc,ac,emin,rmin2,ipon,&s.haspar[0][0],s.potshape,
s.swi_on, s.swi_off, strlen(argv[0]));
  }
  if (argv[1] != NULL) {
    if (single) {
      s.nlig0 = 1;
      int dmmy=0,dmmy2=0;
      read_single_pdb_(
       MAXLIG, TOTMAXRES, TOTMAXATOM, MAXATOM,
       argv[1],s.kai,s.tyi,s.rgi,s.iei,s.x,s.iaci,s.xlai,
       s.icop,s.we,s.we0,s.chai,s.ncop,s.nmaxco,s.natco,
       s.nlig0,s.nres,s.natom,s.n3atom,s.nall,s.nall3,s.ieins,s.ieins3,
       dmmy, dmmy2,
       strlen(argv[1])
      );
    }
    else if (argc == 3) {
      s.nlig0 = 2;
      read_two_pdbs_(
       MAXLIG, TOTMAXRES, TOTMAXATOM, MAXATOM,
       argv[1],argv[2],s.kai,s.tyi,s.rgi,s.iei,s.x,s.iaci,s.xlai,
       s.icop,s.we,s.we0,s.chai,s.ncop,s.nmaxco,s.natco,
       s.nres,s.natom,s.n3atom,s.nall,s.nall3,s.ieins,s.ieins3,
       strlen(argv[1]),strlen(argv[2])
      );
    }
    else {  
      s.nlig0 = 0;
      read_one_pdb_(
       MAXLIG, TOTMAXRES, TOTMAXATOM, MAXATOM,
       argv[1],s.kai,s.tyi,s.rgi,s.iei,s.x,s.iaci,s.xlai,
       s.icop,s.we,s.we0,s.chai,s.ncop,s.nmaxco,s.natco,
       s.nlig0,s.nres,s.natom,s.n3atom,s.nall,s.nall3,s.ieins,s.ieins3,
       strlen(argv[1])
      );
    } 
    s.nlig = s.nlig0;
    memcpy(s.xori0, s.x, 3*TOTMAXATOM*sizeof(double)); 
    for (int n = 0; n < MAXLIG; n++) {
      double *fr = s.forcerotation[n];
      fr[0] = 1;
      fr[4] = 1;
      fr[8] = 1;
    }
  }
    
  for (i=0;i<MAXLIG;i++) {
    s.nhm[i]=0;
    s.nihm[i]=0;
    s.nr_symcopies[i]=1;
    s.symcopies[i][0] = i;
  }  
  memset(s.grids, 0, MAXLIG*sizeof(Grid *));
  return cartstatehandle;
}

CartState &cartstate_get(int handle) {
  return *cartstates[handle-9990];
}  

inline void _cartstate_pivotize(CartState &cartstate) {
  for (int n=0; n < cartstate.nlig; n++) {
    int start = 0;
    if (n) start = cartstate.ieins[n-1];
    for (int nn=start; nn < cartstate.ieins[n]; nn++) {
      cartstate.xb[3*nn] = cartstate.x[3*nn] - cartstate.pivot[0][n];
      cartstate.xb[3*nn+1] = cartstate.x[3*nn+1] - cartstate.pivot[1][n];
      cartstate.xb[3*nn+2] = cartstate.x[3*nn+2] - cartstate.pivot[2][n];      
    }
  }
}

extern "C" void cartstate_set_pivot_(const int &handle, double (&pivot)[MAXLIG][3]) {
  CartState &cartstate = *cartstates[handle-9990];
  memcpy(cartstate.pivot, pivot, 3*MAXLIG*sizeof(double));
  _cartstate_pivotize(cartstate);
}

extern "C" void cartstate_pivot_auto_(const int &handle) {
  CartState &cartstate = *cartstates[handle-9990];
  for (int n=0; n < cartstate.nlig; n++) {
    int start = 0;
    if (n) start = cartstate.ieins[n-1];
    double px=0,py=0,pz=0;
    for (int nn=start; nn < cartstate.ieins[n]; nn++) {
      px += cartstate.x[3*nn];
      py += cartstate.x[3*nn+1];
      pz += cartstate.x[3*nn+2];
    }
    int size = cartstate.ieins[n]-start;
    cartstate.pivot[0][n] = px/size;    
    cartstate.pivot[1][n] = py/size;    
    cartstate.pivot[2][n] = pz/size; 
    //printf("pivot %d %.3f %.3f %.3f\n", n, cartstate.pivot[0][n],cartstate.pivot[1][n],cartstate.pivot[2][n]);       
  }
  _cartstate_pivotize(cartstate);
}

extern "C" void cartstate_get_nlig_nhm_(const int &handle,int &nlig, int *&nhm, int *&nihm) {
  CartState &cartstate = *cartstates[handle-9990];
 
  nhm = &(cartstate.nhm[0]);
  nihm = &(cartstate.nihm[0]);
  nlig = cartstate.nlig;
}

extern "C" void cartstate_get_forcerot_(const int &handle,double *&forcefactor, double *&forcerotation) {
  CartState &cartstate = *cartstates[handle-9990];
 
  forcefactor = &(cartstate.forcefactor[0]);
  forcerotation = &(cartstate.forcerotation[0][0]);
  
}

extern "C" void cartstate_get_forces_(const int &handle,double *&f, int &nall3) {
  CartState &cartstate = *cartstates[handle-9990]; 
  nall3 = cartstate.nall3;
  f = &(cartstate.f[0]);
}

extern "C" void cartstate_get_parameters_(const int &handle,double *&rbc, 
double *&rc,double *&ac,double *&emin,double *&rmin2,int *&ipon, int &potshape, 
int &cdie, double &epsilon, float &swi_on, float &swi_off, int *&use_softcore, double *&softcore) 
{
  CartState &cartstate = *cartstates[handle-9990]; 
  rbc = &(cartstate.rbc[0][0]);
  rc = &(cartstate.rc[0][0]);
  ac = &(cartstate.ac[0][0]);
  emin = &(cartstate.emin[0][0]);
  rmin2 = &(cartstate.rmin2[0][0]);
  ipon = &(cartstate.ipon[0][0]);
  potshape = cartstate.potshape;
  cdie = cartstate.cdie;
  epsilon = cartstate.epsilon;
  swi_on = cartstate.swi_on;
  swi_off = cartstate.swi_off;  
  softcore = &cartstate.softcore;
  use_softcore = &cartstate.use_softcore;
}

extern "C" void cartstate_get_pivot_(const int &handle,double *&pivot) {
  CartState &cartstate = *cartstates[handle-9990]; 
  pivot = &(cartstate.pivot[0][0]);
}

extern "C" void cartstate_get_nrens_(const int &handle,int *&nrens) {
  CartState &cartstate = *cartstates[handle-9990]; 
  nrens = cartstate.nrens;
}

extern "C" void cartstate_get_morphing_(const int &handle,int *&morphing) {
  CartState &cartstate = *cartstates[handle-9990]; 
  morphing = cartstate.morphing;
}

extern "C" void cartstate_get_has_locrests_(const int &handle,int *&has_locrests) {
  CartState &cartstate = *cartstates[handle-9990]; 
  has_locrests = cartstate.has_locrests;
}

extern "C" void cartstate_get_ensd_(const int &handle,
  const int &ligand,
  const int &ens,
  double *&ensd,
  const double &morph,
  double &cmorph,
  double *&cmorphd  
  ) 
  {
  CartState &cartstate = *cartstates[handle-9990]; 
  ensd = cartstate.ensd[ligand][ens-1];
  cmorph = -1;
  int ccmorph = 0;
  if (morph >= 0 && morph < cartstate.nrens[ligand]) {
    ccmorph = int(morph);
    cmorph = morph - ccmorph;
    ensd = cartstate.ensd[ligand][ccmorph];
    cmorphd = cartstate.morphd[ligand][ccmorph];
  }
}


extern "C" void cartstate_f_write_pdb_(
  const int &handle,
  int &nlig, int *&kai, char4 *&tyi, char4 *&rgi, int *&iei, double *&x,
  int *&iaci, double *&xlai, int *&icop, double *&we, int *&ieins) 
{
  CartState &cartstate = *cartstates[handle-9990];

  nlig = cartstate.nlig;
  kai = &(cartstate.kai[0]);
  tyi = &(cartstate.tyi[0]);
  rgi = &(cartstate.rgi[0]);
  iei = &(cartstate.iei[0]);
  x   = &(cartstate.x[0]);
  iaci = &(cartstate.iaci[0]);
  xlai = &(cartstate.xlai[0]);
  icop = &(cartstate.icop[0]);
  we = &(cartstate.we[0]);
  ieins = &(cartstate.ieins[0]);
}

extern "C" void cartstate_f_rotdeform_(
  const int &handle,
  int *&nhm, int *&nihm, int *&ieins, double *&eig, int *&index_eig, double *& index_val,double *&pivot,
  double *&xb, double *&x,
  double *&xori, double *&xori0) 
{

  CartState &cartstate = *cartstates[handle-9990];

  nhm = &(cartstate.nhm[0]);
  nihm = &(cartstate.nihm[0]);
  ieins = &(cartstate.ieins[0]);
  eig = &(cartstate.eig[0][0][0]);
  index_eig = &(cartstate.index_eig[0][0][0]);
  index_val = &(cartstate.index_val[0][0][0]);
  pivot = &(cartstate.pivot[0][0]);
  xb   = &(cartstate.xb[0]);
  x   = &(cartstate.x[0]);
  xori   = &(cartstate.xori[0]);
  xori0   = &(cartstate.xori0[0]);    
}

extern "C" void cartstate_f_pairenergy_(const int &handle,
  int *&nhm, int *&nihm, int *&ieins, double *&eig,
  int *&index_eig, double *&index_val,
  double *&xb, double *&x,double *&xori, double *&xori0) 
{   
  CartState &cartstate = *cartstates[handle-9990];

  nhm = &(cartstate.nhm[0]);
  nihm = &(cartstate.nihm[0]);
  ieins = &(cartstate.ieins[0]);
  eig = &(cartstate.eig[0][0][0]);
  index_eig = &(cartstate.index_eig[0][0][0]);
  index_val = &(cartstate.index_val[0][0][0]);
  xb   = &(cartstate.xb[0]);
  x   = &(cartstate.x[0]);  
  xori   = &(cartstate.xori[0]);
  xori0   = &(cartstate.xori0[0]);  
}

extern "C" void cartstate_f_globalenergy_(const int &handle,
  int &nlig, int &nall, int &nall3, 
  double &morph_fconstant, int *&nhm, int *&nihm, int *&ieins, double *&eig, double *&val,
  int *&index_eig, double *&index_val,double *&xb, double *&x,double *&xori, double *&xori0,
  double *&f, double *&pivot, int *&natom, int *&iaci_old) 
{   
  CartState &cartstate = *cartstates[handle-9990];
  nlig = cartstate.nlig;
  nall = cartstate.nall;
  nall3 = cartstate.nall3;
  morph_fconstant = cartstate.morph_fconstant;
  nhm = &(cartstate.nhm[0]);
  nihm = &(cartstate.nihm[0]);
  ieins = &(cartstate.ieins[0]);
  eig = &(cartstate.eig[0][0][0]);  
  val = &(cartstate.val[0][0]);
  index_eig =  &(cartstate.index_eig[0][0][0]);
  index_val = &(cartstate.index_val[0][0][0]);
  xb   = &(cartstate.xb[0]);
  x   = &(cartstate.x[0]);  
  xori   = &(cartstate.xori[0]);
  xori0   = &(cartstate.xori0[0]);  
  f = &(cartstate.f[0]);
  pivot = &(cartstate.pivot[0][0]);
  natom = &(cartstate.natom[0]);
  iaci_old = &(cartstate.iaci_old[0]);
}

extern "C" void cartstate_f_disre_(
  const int &handle,
  int &nlig, double *&pivot) 
{
  CartState &cartstate = *cartstates[handle-9990];

  nlig = cartstate.nlig;
  pivot = &(cartstate.pivot[0][0]);
}

extern "C" void cartstate_select_ligand_(
 const int &handle,
 const int &ligand,
 int &natom,
 int &nres,
 int *&iei,
 double *&x,
 double *&f,
 double *pivot,
 int *&iaci,
 int *&icop,
 double *&we,
 double *&chai,
 int *&ncop,
 int *&nmaxco,
 int *&natco,
 const int &ori
 ) {
  CartState &cartstate = *cartstates[handle-9990];

  natom = cartstate.natom[ligand];
  nres = cartstate.nres[ligand];
  if (ligand > 0) nres -= cartstate.nres[ligand-1];
 
  int start = 0;
  if (ligand > 0) start = cartstate.ieins[ligand-1];
  iei = &(cartstate.iei[start]);
  if (ori) {
    x = &(cartstate.xori[3*start]);
  }
  else {
    x = &(cartstate.x[3*start]);
  }
  f = &(cartstate.f[3*start]);
  pivot[0] = cartstate.pivot[0][ligand];
  pivot[1] = cartstate.pivot[1][ligand];
  pivot[2] = cartstate.pivot[2][ligand];
  iaci = &(cartstate.iaci[start]);
  icop = &(cartstate.icop[start]);
  we = &(cartstate.we[start]);
  chai = &(cartstate.chai[start]);
  
  int rstart = 0;
  if (ligand > 0) rstart = cartstate.nres[ligand-1];
  nmaxco = &(cartstate.nmaxco[rstart]);
  natco = &(cartstate.natco[rstart]);
  ncop = &(cartstate.ncop[rstart][0][0]);

} 

extern "C" void cartstate_select_ligand2_(
 const int &handle,
 const int &ligand,
 double *&x,
 double *&f
 ) {
  CartState &cartstate = *cartstates[handle-9990]; 
  int start = 0;
  if (ligand > 0) start = cartstate.ieins[ligand-1];
  x = &(cartstate.x[3*start]);
  f = &(cartstate.f[3*start]);
} 

extern "C" void cartstate_get_haspar_(const int &handle,iParameters *&haspar) {
  CartState &cartstate = *cartstates[handle-9990]; 
  haspar = &cartstate.haspar;
}  

extern "C" void cartstate_translate_atomtypes_(const int &handle) {
  CartState &cartstate = *cartstates[handle-9990]; 
  int *transtable = cartstate.transtable;
  int transtable0[MAXATOMTYPES];
  if (transtable == NULL) {
    bool has_32 = 0;
	 transtable = &transtable0[0];
    for (int n = 0; n < MAXATOMTYPES; n++) transtable[n]=n+1;
    transtable[MAXATOMTYPES-1] = 0;
    for (int j = 0; j < MAXATOMTYPES; j++) {
      if (cartstate.haspar[32-1][j]) {
        has_32 = 1;
	break;
      }
    }
    if (!has_32) transtable[32-1] = 0;
  }

  /*
  for (int n = 85; n < 98; n++){
    for (int nn = 0; nn < 98; nn++){
      printf("%d %d %.3f %.3f\n", n+1,nn+1,cartstate.ac[n][nn], cartstate.rc[n][nn]);
    }
  }
  exit(0);
  */

  //make the translation
  for (int n = 0; n < cartstate.nall; n++) {
    cartstate.iaci_old[n] = cartstate.iaci[n];
    int res = cartstate.iaci[n];
    if (res <= 0 || res > MAXATOMTYPES) {
      fprintf(stderr, "Atom %d has invalid type %d\n", n+1, res);
      exit(1);
    }
    cartstate.iaci[n] = transtable[res-1];
  }  
}
extern "C" void cartstate_apply_lambda_(const int  &cartstatehandle) {
  CartState &cartstate = cartstate_get(cartstatehandle);
  if (cartstate.use_lambda) {
    for (int n = 0; n < cartstate.nall; n++) {
      int atomtype = cartstate.iaci[n];
      double f = cartstate.lambda;
      if (atomtype >= 90) f = 1 - cartstate.lambda;
      cartstate.we[n] = cartstate.we0[n] * f;
    }
  }
}


extern "C" void cartstate_apply_epsilon_(const int  &cartstatehandle) {
  CartState &cartstate = cartstate_get(cartstatehandle);
  apply_permi_(TOTMAXATOM, cartstate.nall, cartstate.chai, cartstate.epsilon);
}
   
bool exists(const char *f) {
  FILE *fil = fopen(f, "r");
  if (fil == NULL) return 0;
  else {
    fclose(fil);
    return 1;
  }
}

extern "C" void cartstate_prepare_axsym_(const int &cartstatehandle) {
  CartState &c = cartstate_get(cartstatehandle);
  prepare_axsym_cartstate(c);
}
extern "C" void apply_axsym_(
 const int &cartstatehandle, 
 double *morph, 
 int *ens, 
 double *phi,
 double *ssi,
 double *rot,
 double *xa,
 double *ya,
 double *za,
 modes2 &dlig, 
 coors2 &locrests
)
{
  CartState &c = cartstate_get(cartstatehandle);
  apply_axsym(
   c.nr_symtrans,
   c.symtrans,

   morph, 
   ens, 
   phi,
   ssi,
   rot,
   xa,
   ya,
   za,
   dlig, 
   c.has_locrests,
   locrests
  ); 
}
