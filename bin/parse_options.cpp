#include "state.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <stdexcept>

extern bool exists(const char *);
extern void parse_restraintfile(MiniState &ms, const char *restfile);
extern "C" void define_atomdensitygrid_(float voxelsize, int dimension, float forceconstant);
extern "C" void ignoreweights_atomdensitygrid_();
extern "C" void define_atomdensitymask_(char *maskfile, float forceconstant);

extern void read_ens(int cartstatehandle, int ligand, char *ensfile, bool strict, bool morphing);

extern CartState &cartstate_get(int handle);
extern MiniState &ministate_get(int handle);

void lambda_usage() {
 fprintf(stderr, "--lambda option usage: --lambda <value between 0 (ultra-coarse-grained) and 1 (coarse-grained)\n");
  exit(1);
}

void grid_usage() {
 fprintf(stderr, "--grid/--torquegrid option usage: --grid/--torquegrid <ligand nr> <file name/ligand number>\n");
  exit(1);
}

void ens_usage() {
 fprintf(stderr, "--ens option usage: --ens <ligand nr> <file name>\n");
  exit(1);
}

void morph_usage() {
 fprintf(stderr, "--morph option usage: --morph <ligand nr> <ensemble file name>\n");
  exit(1);
}

void sym_usage() {
 fprintf(stderr, "--sym option usage: --sym <symmetry type> <ligand number> <ligand number> ...\n");
  exit(1); 
}

void axsym_usage() {
  fprintf(stderr, "--axsym option usage: --axsym <ligand> <symmetry>\n  <axis x> <axis y> <axis z>\n  <origin x> <origin y> <origin z>\n");
  exit(1);
}

void ncsym_usage() {
  fprintf(stderr, "--ncsym option usage: --ncsym <ligand> <angle in degrees>\n  <axis x> <axis y> <axis z>\n  <origin x> <origin y> <origin z>\n");
  exit(1);
}

void mctemp_usage() {
 fprintf(stderr, "--mctemp option usage: --mctemp <temperature in KT>\n");
  exit(1);
}

void mcmtemp_usage() {
 fprintf(stderr, "--mcmtemp option usage: --mcmtemp <maximal energy difference>\n");
  exit(1);
}


void epsilon_usage() {
 fprintf(stderr, "--epsilon option usage: --epsilon <dielectric constant; 1=vacuum>\n");
  exit(1);
}

void morph_fconstant_usage() {
 fprintf(stderr, "--morph-fconstant option usage: --morph-fconstant <morphing force constant>\n");
  exit(1);
}

void mcensprob_usage() {
 fprintf(stderr, "--mcensprob option usage: --mcensprob <probability>\n");
  exit(1);
}

void mcscalerot_usage() {
 fprintf(stderr, "--mcscalerot option usage: --mcscalerot <step size in radians>\n");
  exit(1);
}

void mcscalecenter_usage() {
 fprintf(stderr, "--mcscalecenter option usage: --mcscalecenter <step size in A>\n");
  exit(1);
}

void mcscalemode_usage() {
 fprintf(stderr, "--mcscalemode option usage: --mcscalemode <step size in mode A>\n");
  exit(1);
}

void rcut_usage() {
 fprintf(stderr, "--rcut option usage: --rcut <distance squared>\n");
  exit(1);
}

void gravity_usage() {
 fprintf(stderr, "--gravity option usage: --gravity <1,2,3,4 or 5>\n");
  exit(1);
}

void locrest_usage() {
 fprintf(stderr, "--locrest option usage: --locrest <ligand number>\n");
  exit(1);
}

void restweight_usage() {
 fprintf(stderr, "--restweight option usage: --restweight <restraints weight>\n");
  exit(1);
}

void rstk_usage() {
 fprintf(stderr, "--rstk option usage: --rstk <gravity constant>\n");
  exit(1);
}

void vmax_usage() {
 fprintf(stderr, "--vmax option usage: --vmax <maximum number of steps>\n");
  exit(1);
}

void mcmax_usage() {
 fprintf(stderr, "--mcmax option usage: --mcmax <maximum number of MC steps>\n");
  exit(1);
}

void rest_usage() {
 fprintf(stderr, "--rest option usage: --rest <restraint file>\n");
  exit(1);
}

void modes_usage() {
 fprintf(stderr, "--modes option usage: --modes <harmonic modes file>\n");
  exit(1);
}

void imodes_usage() {
	fprintf(stderr, "--imodes option usage: --imodes <index modes file>\n");
	fprintf(stderr, "Use for flexible interface residues/ loops\n");
	  exit(1);

}

void atomdensitygrid_usage() {
 fprintf(stderr, "--atomdensitygrid option usage: --atomdensitygrid <voxel size> <dimension> <force constant>\n");
  exit(1);
}


void atomdensitymask_usage() {
 fprintf(stderr, "--atomdensitymask option usage: --atomdensitymask <mask file> <force constant>\n");
  exit(1);
}

void parse_options(int ministatehandle, int cartstatehandle, int nlig, int argc, char *argv[]) {
  MiniState &ms = ministate_get(ministatehandle);
  CartState &c = cartstate_get(cartstatehandle);
  for (int n = 0; n < argc; n++) {
    char *arg = argv[n];
    if (!strcmp(arg,"--mc")) {
      ms.imc = 1;
    }
    else if (!strcmp(arg, "--mcm")) {
    	ms.imc = 2;
    }
    else if (!strcmp(arg,"--lambda")) {
      if (argc-n < 2) lambda_usage();    
      double lambda = atof(argv[n+1]);
      if (lambda < 0 || lambda > 1) lambda_usage();
      c.use_lambda = 1;
      c.lambda = lambda;
      n += 1;
    }        
    else if (!strcmp(arg,"--morph-fconstant")) {
      if (argc-n < 2) morph_fconstant_usage();    
      double morph_fconstant = atof(argv[n+1]);
      if (morph_fconstant < 0) morph_fconstant_usage();
      c.morph_fconstant = morph_fconstant;
      n += 1;
    }        
    else if (!strcmp(arg,"--epsilon")) {
      if (argc-n < 2) epsilon_usage();    
      double epsilon = atof(argv[n+1]);
      if (epsilon <= 0) epsilon_usage();
      c.epsilon = epsilon;
      n += 1;
    }    
    else if (!strcmp(arg,"--mctemp")) {
      if (argc-n < 2) mctemp_usage();    
      double mctemp = atof(argv[n+1]);
      if (mctemp < 0) mctemp_usage();
      ms.mctemp = mctemp;
      n += 1;
    }    
    else if (!strcmp(arg,"--mcmtemp")) {
      if (argc-n < 2) mcmtemp_usage();
      double mcmtemp = atof(argv[n+1]);
      if (mcmtemp < 0) mcmtemp_usage();
      ms.mcmtemp = mcmtemp;
      n += 1;
    }
    else if (!strcmp(arg,"--mcensprob")) {
      if (argc-n < 2) mcensprob_usage();    
      double mcensprob = atof(argv[n+1]);
      if (mcensprob < 0 || mcensprob > 1) mcensprob_usage();
      ms.mcensprob = mcensprob;
      n += 1;
    }    
    else if (!strcmp(arg,"--mcscalerot")) {
      if (argc-n < 2) mcscalerot_usage();    
      double mcscalerot = atof(argv[n+1]);
      if (mcscalerot < 0) mcscalerot_usage();
      ms.mcscalerot = mcscalerot;
      n += 1;
    }    
    else if (!strcmp(arg,"--mcscalecenter")) {
      if (argc-n < 2) mcscalecenter_usage();    
      double mcscalecenter = atof(argv[n+1]);
      if (mcscalecenter < 0) mcscalecenter_usage();
      ms.mcscalecenter = mcscalecenter;
      n += 1;
    }    
    else if (!strcmp(arg,"--mcscalemode")) {
      if (argc-n < 2) mcscalemode_usage();    
      double mcscalemode = atof(argv[n+1]);
      if (mcscalemode < 0) mcscalemode_usage();
      ms.mcscalemode = mcscalemode;
      n += 1;
    }    
    else if (!strcmp(arg,"--score")) {
      ms.iscore = 1;
    }
    else if (!strcmp(arg,"--traj")) {
      ms.iscore = 2;
    }
    else if (!strcmp(arg,"--cdie")) {
      c.cdie = 1;
    }
    else if (!strcmp(arg,"--ghost")) {
      ms.ghost = 1;
    }
    else if (!strcmp(arg,"--fix-receptor")) {
      ms.fixre = 1;
    }
    else if (!strcmp(arg,"--ghost-ligands")) {
      ms.ghost_ligands = 1;
    }    
    else if (!strcmp(arg,"--only-rot")) {
      ms.itra = 0;
      ms.iori = 1;
    }
    else if (!strcmp(arg,"--only-trans")) {
      ms.itra = 1;
      ms.iori = 0;
    }    
    else if (!strcmp(arg,"--only-flex")) {
      ms.itra = 0;
      ms.iori = 0;
    }    
    else if (!strcmp(arg,"--use-softcore")) {

	c.use_softcore = 1; //MANU softcore wird benutzt
	n += 1;
	try{
		c.softcore = atof(argv[n]); //gewünschten softcore aus commandline ziehen und in cartstate Variable softcore schreiben	
	}
	catch(const std::invalid_argument &e){
		n = n-1;
		c.softcore = 5.0;	//Standartwert 5 in Cartstatevariable Softcore schreiben	
	}   
      }   
    else if ((!strcmp(arg,"--grid"))||(!strcmp(arg,"--torquegrid"))) {
      bool torquegrid = 0;
      if (!strcmp(arg,"--torquegrid")) torquegrid = 1;
      if (argc-n < 3) grid_usage();
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) grid_usage();
      char *gridf = argv[n+2];
      if (!exists(gridf)) {
        char *endptr;
        int lig_old = strtol(gridf,&endptr, 0);
        if (lig_old == 0 || endptr-gridf < (unsigned int) strlen(gridf)) {        
          fprintf(stderr, "Grid file %s does not exist\n", gridf);
  	  grid_usage();          
        }
        if (lig_old < 0 || lig_old > nlig || c.grids[lig_old-1] == NULL) {
          fprintf(stderr, "Grid index %d does not point to an existing grid\n", lig_old);
  	  grid_usage();          
        }
        c.grids[lig-1] = c.grids[lig_old-1];
        n += 2;
        continue;
      }
      Grid *g = new Grid;
      if (torquegrid) {
        g->read_torque(gridf);
      }  
      else {  
        g->read_std(gridf);
      }  
      if (g->natoms != c.natom[lig-1]) {
        fprintf(stderr, "Wrong number of atoms for ligand %d:\n  Grid file %s: %d, PDB file: %d\n",lig,gridf,g->natoms,c.natom[lig-1]);
	exit(1);
      }
      c.grids[lig-1] = g;
      n += 2;
    }
    else if (!strcmp(arg,"--sym")) {
      if (argc-n < 4) sym_usage();
      int type = atoi(argv[n+1]);      
      int nrsym = type;
      if (type == -2) {
        nrsym = 4;      
      }      
      else {      
        if (type < 2 || type > MAXLIG) sym_usage();
      }
      //printf("%d %d %d\n", type, nrsym, MAXLIG);
      if (argc-n < nrsym + 2) sym_usage();
      if (c.nsym >= MAXLIG-1) {
        fprintf(stderr, "Too many symmetries\n");
        exit(1);
      }
      c.symtypes[c.nsym] = type;
      for (int nn = 0; nn  < nrsym; nn++) {
        int lig = atoi(argv[n+2+nn]);
        if (lig < 1 || lig > nlig) sym_usage();
        c.sym[c.nsym][nn] = lig;
      }
      c.nsym++;
      ms.has_globalenergy = 1;
      n += 1+nrsym;
    }
    else if (!strcmp(arg,"--axsym")) {
      if (argc-n < 3) axsym_usage();
      AxSymmetry &sym = c.axsyms[c.nr_axsyms];  
      sym.ligand = atoi(argv[n+1]);
      sym.symtype = atoi(argv[n+2]);
      sym.angle = 0;
      sym.axis[0] = atof(argv[n+3]);
      sym.axis[1] = atof(argv[n+4]);
      sym.axis[2] = atof(argv[n+5]);
      sym.origin[0] = atof(argv[n+6]);
      sym.origin[1] = atof(argv[n+7]);
      sym.origin[2] = atof(argv[n+8]);  
      c.nr_axsyms++;
      n += 8;
    }
    else if (!strcmp(arg,"--ncsym")) {
      if (argc-n < 3) ncsym_usage();
      AxSymmetry &sym = c.axsyms[c.nr_axsyms];  
      sym.ligand = atoi(argv[n+1]);
      sym.symtype = 0;
      sym.angle = atof(argv[n+2]);
      sym.axis[0] = atof(argv[n+3]);
      sym.axis[1] = atof(argv[n+4]);
      sym.axis[2] = atof(argv[n+5]);
      sym.origin[0] = atof(argv[n+6]);
      sym.origin[1] = atof(argv[n+7]);
      sym.origin[2] = atof(argv[n+8]);  
      c.nr_axsyms++;
      n += 8;
    }    
    else if (!strcmp(arg,"--ens")) {
      if (argc-n < 3) ens_usage();
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) ens_usage();
      char *ensf = argv[n+2];
      if (!exists(ensf)) {
        fprintf(stderr, "Ensemble file %s does not exist\n", ensf);
	ens_usage();
      }
      if (c.nrens[lig-1]) {
        fprintf(stderr, "Ligand %d can have only one ensemble/morphing specification\n", lig);
        exit(1);
      }
      read_ens(cartstatehandle,lig-1,ensf,1,0);
      n += 2;
    }
    else if (!strcmp(arg,"--morph")) {
      if (argc-n < 3) morph_usage();
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) morph_usage();
      char *ensf = argv[n+2];
      if (!exists(ensf)) {
        fprintf(stderr, "Ensemble file %s does not exist\n", ensf);
	morph_usage();
      }
      if (c.nrens[lig-1]) {
        fprintf(stderr, "Ligand %d can have only one ensemble/morphing specification\n", lig);
        exit(1);
      }
      c.morphing[lig-1] = 1;
      read_ens(cartstatehandle,lig-1,ensf,1,1);
      ms.has_globalenergy = 1;
      n += 2;
    }
    else if (!strcmp(arg,"--rcut")) {
      if (argc-n < 2) rcut_usage();    
      double rcut = atof(argv[n+1]);
      if (rcut <= 0) rcut_usage();
      ms.rcut = rcut;
      n += 1;
    }
    else if (!strcmp(arg,"--gravity")) {
      if (argc-n < 2) gravity_usage();    
      int gravity = atoi(argv[n+1]);
      if (gravity <= 0 || gravity > 5) gravity_usage();
      ms.gravity = gravity;
      n += 1;
    }
    else if (!strcmp(arg,"--locrest")) {
      if (argc-n < 2) locrest_usage();    
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) locrest_usage();
      c.has_locrests[lig-1] = 1;
      n += 1;
    }
    else if (!strcmp(arg,"--rstk")) {
      if (argc-n < 2) rstk_usage();    
      double rstk = atof(argv[n+1]);
      ms.rstk = rstk;
      n += 1;
    }    
    else if (!strcmp(arg,"--restweight")) {
      if (argc-n < 2) restweight_usage();    
      double restweight = atof(argv[n+1]);
      ms.restweight = restweight;
      n += 1;
    }        
    else if (!strcmp(arg,"--vmax")) {
      if (argc-n < 2) vmax_usage();    
      int vmax = atoi(argv[n+1]);
      if (vmax <= 0) vmax_usage();
      ms.ivmax = vmax;
      n += 1;
    }
    else if (!strcmp(arg,"--mcmax")) {
      if (argc-n < 2) mcmax_usage();
      int mcmax = atoi(argv[n+1]);
      if (mcmax <= 0) mcmax_usage();
      ms.imcmax = mcmax;
      n += 1;
    }
    
    else if (!strcmp(arg,"--rest")) {
      if (argc-n < 2) rest_usage();    
      char *restf = argv[n+1];
      if (!exists(restf)) {
        fprintf(stderr, "Restraint file %s does not exist\n", restf);
	rest_usage();
      }
      parse_restraintfile(ms, restf);
      ms.has_globalenergy = 1;
      n += 1;
    }
    else if (!strcmp(arg,"--modes")) {
      if (argc-n < 2) modes_usage();    
      char *hmf = argv[n+1];
      if (!exists(hmf)) {
        fprintf(stderr, "Modes file %s does not exist\n", hmf);
	modes_usage();
      }
      const int multi = 1;
      read_hm_(hmf, "ligand", c.nlig, c.natom, c.nhm, c.val, (double *) c.eig, multi, strlen(hmf), strlen("ligand"));      
      n += 1;
      ms.ieig = 1;
      if (ms.has_globalenergy == 0) ms.has_globalenergy = 2;
    }
    else if (!strcmp(arg,"--imodes")) {
      if (argc-n < 2) modes_usage();
      char *hmf = argv[n+1];
      if (!exists(hmf)) {
        fprintf(stderr, "Index modes file %s does not exist\n", hmf);
	imodes_usage();
      }
      const int multi = 1;
      read_indexmode_(hmf, "ligand", c.nlig, c.nihm, (int *) c.index_eig, (double *) c.index_val, multi, strlen(hmf), strlen("ligand"));
      n += 1;
      ms.iindex = 1;
    }
    else if (!strcmp(arg,"--atomdensitygrid")) {
      if (argc-n < 4) atomdensitygrid_usage();    
      float voxelsize = atof(argv[n+1]);
      if (voxelsize <= 0) atomdensitygrid_usage();    
      int dimension = atoi(argv[n+2]);
      if (dimension <= 0 || dimension >= 4000) atomdensitygrid_usage();    
      float forceconstant = atof(argv[n+3]);
      if (forceconstant <= 0) atomdensitygrid_usage();
      define_atomdensitygrid_(voxelsize, dimension, forceconstant);
      ms.has_globalenergy = 1;
      n += 3;
    } 
    else if (!strcmp(arg,"--ignoreatomdensityweights")) {
      ignoreweights_atomdensitygrid_();
    }
    else if (!strcmp(arg,"--atomdensitymask")) {
      if (argc-n < 3) atomdensitymask_usage();    
      char *maskfile = argv[n+1];
      if (!exists(maskfile)) {
        fprintf(stderr, "Mask file %s does not exist\n", maskfile);
        exit(1);
      }      
      float forceconstant = atof(argv[n+2]);
      if (forceconstant <= 0) atomdensitymask_usage();
      define_atomdensitymask_(maskfile, forceconstant);
      ms.has_globalenergy = 1;
      n += 2;
    }        
    else {
      fprintf(stderr, "Unknown option %s\n", arg);
      exit(1);
    }
  }
}
