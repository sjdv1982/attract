#include "axsym.h"
#include "state.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

extern CartState &cartstate_get(int handle);

extern "C" void ministate_check_parameters_(const int &ministatehandle, const int  &cartstatehandle);

MiniState *ministates[100]; 
int ministatesize = 0;

extern "C" int ministate_new_() {
  MiniState *ms0 = new MiniState[1];
  MiniState &ms = *ms0;
  memset(&ms, 0, sizeof(MiniState));
  ministates[ministatesize] = &ms;
  ministatesize++;
  //default settings
  ms.imc = 0;
  ms.mctemp = 1.0;
  ms.mcmtemp = 15.0;
  ms.mcscalerot = 0.05;
  ms.mcscalecenter = 0.1;
  ms.mcscalemode = 3;
  ms.mcensprob = 0.05;
  ms.iscore = 0;
  ms.ivmax = 100;
  ms.imcmax = 100;
  ms.iori = 1;
  ms.itra = 1; 
  ms.ieig = 0;
  ms.iindex = 0;
  ms.irst = 0; 
  ms.fixre = 0; 
  ms.rcut = 1500;
  ms.nr_restraints = 0;   
  ms.has_globalenergy = 0;
  ms.gravity = 0;
  ms.rstk = 0.2; //for harmonic restraints...
  ms.ghost = 0;
  ms.ghost_ligands = 0;
  return ministatesize-1; 

}

MiniState &ministate_get(int handle) {
  return *ministates[handle];
}

extern "C" void ministate_ghost_(const int &handle, int &ghost) {
  MiniState &ms = *ministates[handle];
  ghost = ms.ghost;
}

extern "C" void ministate_iscore_imc_(const int &handle, int &iscore, int &imc) {
  MiniState &ms = *ministates[handle];
  iscore = ms.iscore;
  imc = ms.imc;
}

extern "C" void ministate_f_minfor_(const int &handle, int &iscore, int &ivmax, int &iori, int &itra, int &ieig, int &iindex, int &fixre, int &gridmode) {
  MiniState &ms = *ministates[handle];
  iscore=ms.iscore; ivmax = ms.ivmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig;iindex = ms.iindex;
  fixre = ms.fixre; gridmode = ms.gridmode;
}

extern "C" void ministate_f_minfor_min_(const int &handle, int &iscore, int &ivmax, int &iori, int &itra, int &ieig, int &iindex, int &fixre, int &gridmode) {
  MiniState &ms = *ministates[handle];
  iscore=ms.iscore; ivmax = ms.ivmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig;iindex = ms.iindex;
  fixre = ms.fixre; gridmode = ms.gridmode;
}

extern "C" void ministate_f_monte_(const int &handle, int &iscore, int &imcmax, int &iori, int &itra, int &ieig, int &iindex, int &fixre, int &gridmode, double &mctemp, double &mcscalerot, double &mcscalecenter, double &mcscalemode, double &mcensprob) {
  MiniState &ms = *ministates[handle];
  iscore=ms.iscore; imcmax = ms.imcmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig; iindex=ms.iindex;
  fixre = ms.fixre; gridmode = ms.gridmode; 
  mctemp = ms.mctemp; mcscalerot = ms.mcscalerot; mcscalecenter = ms.mcscalecenter; mcscalemode = ms.mcscalemode; mcensprob = ms.mcensprob;
}

extern "C" void ministate_f_monte_min_(const int &handle, int &iscore, int &imcmax, int &iori, int &itra, int &ieig, int &iindex, int &fixre, int &gridmode, double &mctemp, double &mcmtemp, double &mcscalerot, double &mcscalecenter, double &mcscalemode, double &mcensprob) {
  MiniState &ms = *ministates[handle];
  iscore=ms.iscore; imcmax = ms.imcmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig; iindex=ms.iindex;
  fixre = ms.fixre; gridmode = ms.gridmode;
  mctemp = ms.mctemp; mcmtemp = ms.mcmtemp; mcscalerot = ms.mcscalerot; mcscalecenter = ms.mcscalecenter; mcscalemode = ms.mcscalemode; mcensprob = ms.mcensprob;
}

extern "C" void ministate_f_disre_(const int &handle, int &gravity, double &rstk) {
  MiniState &ms = *ministates[handle];
  gravity = ms.gravity;
  rstk = ms.rstk;
}

extern "C" void ministate_has_globalenergy_(const int &handle, int &has_globalenergy) {
  MiniState &ms = *ministates[handle];
  has_globalenergy=ms.has_globalenergy;
}

extern "C" void select_(const int &maxatom, const int &maxres, const int &maxmolpair,
  const int &molpairhandle, const int &cartstatehandle, 
  const double &rcut);

extern "C" void pairgen_(const int &maxatom, const int &maxres, const int &maxmolpair,
  const int &molpairhandle, const int &cartstatehandle, 
  const double &rcut);

extern "C" void ministate_calc_pairlist_(const int &ministatehandle, const int  &cartstatehandle) {
  /*
   * This routine creates molpairs (partner-partner interactions)
   * It checks that all molpairs are grid-accellerated (or none are)
   * It also checks for inconsistencies like receptor modes + torque grids 
   */
  MiniState &ms = *ministates[ministatehandle];
  CartState &cartstate = cartstate_get(cartstatehandle);
  if (ms.pairs != NULL) {
    fprintf(stderr, "Memory leak, ms.pairs is not NULL!");
    exit(1);
  }
  //printf("%d", ms.imc);
  ms.pairs = new MolPair[cartstate.nlig*(cartstate.nlig-1)]; 
  memset(ms.pairs, 0, cartstate.nlig*(cartstate.nlig-1)*sizeof(MolPair));
  int molpairindex = 0;
  /*
   * For every partner-partner interaction, we have to create one or two molpairs
   * The second molpair is only to get the forces on the "receptor" (= the first molecule)
   * Secondary molpairs are only ever created for standard grids (gridmode == 1)
   *  if we are not using grids (or if we are using torque grids), the receptor forces are already taken care of
   *  if we are using Monte Carlo, we don't care about forces and we never create a secondary molpair
   */
  if ((!ms.fixre) && ms.ghost_ligands) {
    fprintf(stderr, "With ghost ligands, you must fix the receptor");
    exit(1);
  }  
  for (int i = 0; i < cartstate.nlig; i++) {
    for (int j = i+1; j < cartstate.nlig; j++) {
      bool two_molpairs = 0;
      
      if (ms.gridmode == 1 && ms.imc == 0) { /*use normal grids, no Monte Carlo*/
	two_molpairs = 1;
        if (ms.fixre) {
          if (i == 0 && (cartstate.nhm[i] == 0 || !ms.ieig) ) {
            two_molpairs = 0;
          }
        }
      }
      else if (ms.gridmode == 2 || (ms.gridmode == 1 && ms.imc) ) { /*use torque grids OR Monte Carlo (we don't care about forces)*/
	two_molpairs = 0;
	/* If you use torque grids, then any partner with modes must be matched by a grid on the other partner */
	if (!ms.imc && ms.ieig && cartstate.nhm[i] && cartstate.nhm[j]) two_molpairs = 1;
      }
      else { /*use no grids at all*/
	two_molpairs = 0;
      }
      
      int m1 = i, m2 = j;      
      int error = 0;
      if (ms.gridmode) {
	
	if (two_molpairs){
	  if (cartstate.grids[i] == NULL || cartstate.grids[j] == NULL) error = 1;
	}
	else {
	  error = 1;
	  if (cartstate.grids[i] != NULL && (ms.imc || !ms.ieig  || cartstate.nhm[i] == 0) ) { 	    
	    error = 0;
	  }  
	  else if (cartstate.grids[j] != NULL && (ms.imc || !ms.ieig || cartstate.nhm[j] == 0) ) {
	    m1 = j;
	    m2 = i;
	    error = 0;	    
	  }
	  if (cartstate.nhm[m1]) error = 2;
	}  	
      }
      if (ms.ghost_ligands && m1 > 0) continue; //in ghost-ligand mode, skip the molpair if the first one isn't the receptor
	
      if (error == 1) {
        printf("ERROR: using a single grid is not possible! Please recheck your input (maybe you are using modes?)\n");
        exit(1);
      }
      else if (error == 2){
        printf("ERROR: when using modes on a ligand the grid for the other ligand has to be supplied!\n");
	exit(1); 
      }
            
      //printf("PAIR %d %d %d %d %d %d\n", i,j, is_grid, ms.gridmode, cartstate.grids[i],cartstate.grids[j]);          
      MolPair &mp = ms.pairs[molpairindex];
      mp.use_energy = 1;
      mp.receptor = m1;
      mp.ligand = m2;      
      if (ms.gridmode) {
        mp.grid = cartstate.grids[m1];
      }
      else {
         mp.pairgen_done = 0;
      }
      molpairindex++;
      
      if (two_molpairs) { 
	MolPair &mp = ms.pairs[molpairindex];
	mp.use_energy = 0;
	mp.receptor = m2;
	mp.ligand = m1;      
        mp.grid = cartstate.grids[m2];
	molpairindex++;	
      }
    }
  }
  ms.npairs = molpairindex;
  ministate_check_parameters_(ministatehandle, cartstatehandle);
}

extern "C" void ministate_free_pairlist_(const int &handle) {
  MiniState &ms = *ministates[handle];
  for (int n = 0; n < ms.npairs; n++) {
    MolPair &mp = ms.pairs[n];
    delete[] mp.iactl;
    delete[] mp.iactr;
    delete[] mp.nonl;
    delete[] mp.nonr;
  }
  delete[] ms.pairs;
  ms.pairs = NULL;
}

extern "C" void ministate_get_molpairhandles_(const int &handle,
 int *mphandles, int &npairs) {
  MiniState &ms = *ministates[handle];
  npairs = ms.npairs;
  for (int n = 0; n < npairs;n++) {
    mphandles[n] = MAXMOLPAIR * (handle) + n + 1;
  }
}
extern "C" void ministate_get_molpairhandle_(const int &ministatehandle,
 const int &receptor, const int &ligand,
 int &molpairhandle) {
  molpairhandle = -1;

  MiniState &ms = *ministates[ministatehandle]; 
  for (int n = 0; n < ms.npairs; n++) {
    MolPair &mp = ms.pairs[n];
    if (mp.receptor == receptor && mp.ligand == ligand) {
      molpairhandle = MAXMOLPAIR * (ministatehandle) + n + 1;
      return;
    }
  }
}

extern "C" void molpair_get_rl_(
 const int &molpairhandle,
 int &receptor,int &ligand,Grid *&gridptr, int &use_energy
) {
  int ministatehandle = int(molpairhandle/MAXMOLPAIR); //rounds down...
  MiniState &ms = *ministates[ministatehandle];
  MolPair &mp = ms.pairs[molpairhandle-MAXMOLPAIR*ministatehandle-1];
  receptor = mp.receptor;
  ligand = mp.ligand;
  gridptr = mp.grid;
  use_energy = mp.use_energy;
}

extern "C" void molpair_get_values_(
 const int &molpairhandle,
 int &receptor,int &ligand,
 int *&iactr,int *&iactl,
 int &nonp, int *&nonr,int *&nonl) {
  int ministatehandle = molpairhandle/MAXMOLPAIR; //rounds down...
  MiniState &ms = *ministates[ministatehandle]; 
  MolPair &mp = ms.pairs[molpairhandle-MAXMOLPAIR*ministatehandle-1];
  receptor = mp.receptor;
  ligand = mp.ligand;
  iactr = &(mp.iactr[0]);
  iactl = &(mp.iactl[0]);
  nonp = mp.nonp;
  nonr = &(mp.nonr[0]);
  nonl = &(mp.nonl[0]);
}

extern "C" void molpair_set_nonp_(
 const int &molpairhandle,
 const int &nonp) {
  int ministatehandle = molpairhandle/MAXMOLPAIR; //rounds down...
  MiniState &ms = *ministates[ministatehandle];
  MolPair &mp = ms.pairs[molpairhandle-MAXMOLPAIR*ministatehandle-1];
  mp.nonp = nonp;
}

extern "C" void molpair_pairgen_(
 const int &molpairhandle,
 const int &cartstatehandle,
 int &nonp
 ) {
  
  //Will generate nonbonded atom pairs, if not already done
 
  int ministatehandle = molpairhandle/MAXMOLPAIR; //rounds down...
  MiniState &ms = *ministates[ministatehandle]; 
  MolPair &mp = ms.pairs[molpairhandle-MAXMOLPAIR*ministatehandle-1];
  if (!mp.pairgen_done) {
    mp.iactl = new int[MAXATOM];
    mp.iactr = new int[MAXATOM];
    mp.nonr = new int[MAXMOLPAIR];
    mp.nonl = new int[MAXMOLPAIR];
    select_(MAXATOM, MAXRES, MAXMOLPAIR, 
      molpairhandle, cartstatehandle, ms.rcut);       
    pairgen_(MAXATOM, MAXRES, MAXMOLPAIR,
      molpairhandle, cartstatehandle, ms.rcut);
    mp.pairgen_done = 1;
  }
  nonp = mp.nonp;
}



extern "C" void ministate_check_parameters_(const int &ministatehandle, const int  &cartstatehandle) {
  MiniState &ms = *ministates[ministatehandle];
  CartState &cartstate = cartstate_get(cartstatehandle);

  //determine if we have all the parameter we need
  //TODO: check the grid alphabets!  
  
  bool used[MAXLIG][MAXATOMTYPES];
  memset(used, 0, sizeof(used));
  for (int n = 0; n < cartstate.nlig; n++) {
    int start = 0;
    if (n > 0) start = cartstate.ieins[n-1];
    for (int nn = start; nn < cartstate.ieins[n]; nn++) {
      int res = cartstate.iaci[nn];
      if (res == 0) continue;
      if (res < 0) {
         fprintf(stderr, "Invalid atom type %d\n", res-1);      
         exit(1);
       }
      used[n][res-1] = 1;
    }
  }
  
  for (int n = 0; n < ms.npairs; n++) {
    MolPair &p = ms.pairs[n];
    //printf("%d %d %d\n",n, p.receptor,p.ligand);
    bool *used1 = &used[p.receptor][0];
    bool *used2 = &used[p.ligand][0];
    for (int i = 0; i < MAXATOMTYPES; i++) {
      if (!used1[i]) continue; 
      for (int ii = 0; ii < MAXATOMTYPES; ii++) {
        if (!used2[ii]) continue;
        if (!cartstate.haspar[i][ii]) {
          fprintf(stderr, 
            "Atom type %d in molecule %d can interact with atom type %d in molecule %d, but no parameters are available\n", i+1,p.receptor+1,ii+1,p.ligand+1);
          exit(1);
        }
      }
    }
  }
}

extern "C" void axsym_fold_grads_(
 const int &ministatehandle,
 const int &cartstatehandle, 
 double *grads0, double *grads, double *morphgrads
) {
  MiniState &ms = *ministates[ministatehandle];
  CartState &cartstate = cartstate_get(cartstatehandle);
  axsym_fold_grads(ms, cartstate, grads0, grads, morphgrads);
}
