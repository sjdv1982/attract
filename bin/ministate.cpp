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
  ministates[ministatesize] = &ms;
  ministatesize++;
  //default settings
  ms.imc = 0;
  ms.mctemp = 3.5;
  ms.mcscalerot = 0.05;
  ms.mcscalecenter = 0.1;
  ms.mcscalemode = 3;
  ms.iscore = 0;
  ms.ivmax = 100; 
  ms.iori = 1;
  ms.itra = 1; 
  ms.ieig = 0;
  ms.irst = 0; 
  ms.fixre = 0; 
  ms.rcut = 1500;
  ms.proxlim = 36;
  ms.proxmax = 200;
  ms.proxmaxtype = 31;
  ms.nr_restraints = 0;
#ifdef TORQUEGRID
  ms.gridmode = 2;
#else 
  ms.gridmode = 1;
#endif   
  ms.has_globalenergy = 0;
  ms.gravity = 0;
  //ms.rstk = 0.015; //for second order restraints...; too much?
  ms.rstk = 0.0015; //for second order restraints...
  return ministatesize-1+7770; 

}

MiniState &ministate_get(int handle) {
  return *ministates[handle-7770];
}

void ministate_iscore_imc_(const int &handle, int &iscore, int &imc) {
  MiniState &ms = *ministates[handle-7770];
  iscore = ms.iscore;
  imc = ms.imc;
}

extern "C" void ministate_f_minfor_(const int &handle, int &iscore, int &ivmax, int &iori, int &itra, int &ieig, int &fixre, int &gridmode) {
  MiniState &ms = *ministates[handle-7770];
  iscore=ms.iscore; ivmax = ms.ivmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig; fixre = ms.fixre; gridmode = ms.gridmode;
}

extern "C" void ministate_f_monte_(const int &handle, int &iscore, int &ivmax, int &iori, int &itra, int &ieig, int &fixre, int &gridmode, double &mctemp, double &mcscalerot, double &mcscalecenter, double &mcscalemode) {
  MiniState &ms = *ministates[handle-7770];
  iscore=ms.iscore; ivmax = ms.ivmax;iori = ms.iori;itra = ms.itra;ieig = ms.ieig; 
  fixre = ms.fixre; gridmode = ms.gridmode; 
  mctemp = ms.mctemp; mcscalerot = ms.mcscalerot; mcscalecenter = ms.mcscalecenter; mcscalemode = ms.mcscalemode;
}

extern "C" void ministate_f_disre_(const int &handle, int &gravity, double &rstk) {
  MiniState &ms = *ministates[handle-7770];
  gravity = ms.gravity;
  rstk = ms.rstk;
}

extern "C" void ministate_has_globalenergy_(const int &handle, int &has_globalenergy) {
  MiniState &ms = *ministates[handle-7770];
  has_globalenergy=ms.has_globalenergy;
}

extern "C" void select_(const int &maxatom, const int &maxres, const int &maxmolpair,
  const int &molpairhandle, const int &cartstatehandle, 
  const double &rcut);

extern "C" void pairgen_(const int &maxatom, const int &maxres, const int &maxmolpair,
  const int &molpairhandle, const int &cartstatehandle, 
  const double &rcut);

extern "C" void ministate_calc_pairlist_(const int &ministatehandle, const int  &cartstatehandle) {
  MiniState &ms = *ministates[ministatehandle-7770];
  CartState &cartstate = cartstate_get(cartstatehandle);
  ms.pairs = new MolPair[cartstate.nlig*(cartstate.nlig-1)]; 
  int molpairindex = 0;
  for (int i0 = 0; i0 < cartstate.nlig; i0++) {
    for (int j0 = 0; j0 < cartstate.nlig; j0++) {
      if (i0 == j0) continue;
      int i = i0, j = j0;
      bool is_grid = 0;
      
      if (ms.gridmode == 1) {
        if (ms.fixre) {
	  if (i == 0) {
	    if (cartstate.grids[i] != NULL) is_grid = 1;
	  }
	  else if (j == 0) {
	    continue;
	  }
	  else if (cartstate.grids[i] != NULL && cartstate.grids[j] != NULL) {
	    is_grid = 1;
	  }
	  else {
	    if (j < i) continue;
	  }	   	  
	}
	else {
	  if (cartstate.grids[i] != NULL && cartstate.grids[j] != NULL) {
	    is_grid = 1;
	  }
	  else {
	    if (j < i) continue;
	  }
	}
      }
      else if (ms.gridmode == 2) {
        if (j < i) continue; 
	//if the receptor has a grid, use it
	if (cartstate.grids[i] != NULL) {
          is_grid = 1;
	}      
	//else if the ligand has a grid, swap receptor and ligand and use it
	else if (cartstate.grids[j] != NULL) {
          i = j0;
	  j = i0;
          is_grid = 1;
	}
      }  
      //printf("PAIR %d %d %d %d %d %d\n", i,j, is_grid, ms.gridmode, cartstate.grids[i],cartstate.grids[j]);          
      MolPair &mp = ms.pairs[molpairindex];
      mp.receptor = i;
      mp.ligand = j;      
      //int molpairhandle = 1000 * ministatehandle + molpairindex + 1;
      if (is_grid) {
        mp.grid = cartstate.grids[i];
      }
      else {
         mp.pairgen_done = 0;
      }
      molpairindex++;
    }
  }
  ms.npairs = molpairindex;
  ministate_check_parameters_(ministatehandle, cartstatehandle);
}

extern "C" void ministate_free_pairlist_(const int &handle) {
  MiniState &ms = *ministates[handle-7770];
  delete[] ms.pairs;
}

extern "C" void ministate_get_molpairhandles_(const int &handle,
 int *mphandles, int &npairs) {
  MiniState &ms = *ministates[handle-7770];
  npairs = ms.npairs;
  for (int n = 0; n < npairs;n++) {
    mphandles[n] = 1000 * handle + n + 1;
  }
}
extern "C" void ministate_get_molpairhandle_(const int &ministatehandle,
 const int &receptor, const int &ligand,
 int &molpairhandle) {
  molpairhandle = -1;

  MiniState &ms = *ministates[ministatehandle-7770]; 
  for (int n = 0; n < ms.npairs; n++) {
    MolPair &mp = ms.pairs[n];
    if (mp.receptor == receptor && mp.ligand == ligand) {
      molpairhandle = 1000 * ministatehandle + n + 1;
      return;
    }
  }
}

extern "C" void molpair_get_rl_(
 const int &molpairhandle,
 int &receptor,int &ligand,Grid *&gridptr
) {
  int ministatehandle = int(molpairhandle/1000); //rounds down...
  MiniState &ms = *ministates[ministatehandle-7770]; 
  MolPair &mp = ms.pairs[molpairhandle-1000*ministatehandle-1];
  receptor = mp.receptor;
  ligand = mp.ligand;
  gridptr = mp.grid;
}

extern "C" void molpair_get_values_(
 const int &molpairhandle,
 int &receptor,int &ligand,
 int *&iactr,int *&iactl,
 int &nonp, int *&nonr,int *&nonl) {
  int ministatehandle = molpairhandle/1000; //rounds down...
  MiniState &ms = *ministates[ministatehandle-7770]; 
  MolPair &mp = ms.pairs[molpairhandle-1000*ministatehandle-1];
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
  int ministatehandle = molpairhandle/1000; //rounds down...
  MiniState &ms = *ministates[ministatehandle-7770]; 
  MolPair &mp = ms.pairs[molpairhandle-1000*ministatehandle-1];
  mp.nonp = nonp;
}

extern "C" void molpair_pairgen_(
 const int &molpairhandle,
 const int &cartstatehandle,
 int &nonp
 ) {
  
  //Will generate nonbonded atom pairs, if not already done
 
  int ministatehandle = molpairhandle/1000; //rounds down...
  MiniState &ms = *ministates[ministatehandle-7770]; 
  MolPair &mp = ms.pairs[molpairhandle-1000*ministatehandle-1];
  if (!mp.pairgen_done) {
    select_(MAXATOM, MAXRES, MAXMOLPAIR, 
      molpairhandle, cartstatehandle, ms.rcut);       
    pairgen_(MAXATOM, MAXRES, MAXMOLPAIR,
      molpairhandle, cartstatehandle, ms.rcut);
    mp.pairgen_done = 1;
  }
  nonp = mp.nonp;
}



extern "C" void ministate_check_parameters_(const int &ministatehandle, const int  &cartstatehandle) {
  MiniState &ms = *ministates[ministatehandle-7770];
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
      //printf("%d\n", res-1);      
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
	    "Atom type %d in molecule %d can interact with atom type %d in molecule %d, but no parameters are available\n", i+1,p.receptor+1,ii+i,p.ligand+1);
	  exit(1);
	}
      }
    }
  }
}
