#include "state.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

extern bool exists(const char *);
extern void parse_restraintfile(MiniState &ms, const char *restfile);
extern "C" void read_densitymaps_(char *densitymapsfile0, int len_densitymapsfile);

extern void read_ens(int cartstatehandle, int ligand, char *ensfile, bool strict);

extern CartState &cartstate_get(int handle);
extern MiniState &ministate_get(int handle);

void grid_usage() {
 fprintf(stderr, "--grid option usage: --grid <ligand nr> <file name>\n");
  exit(1);
}

void ens_usage() {
 fprintf(stderr, "--ens option usage: --ens <ligand nr> <file name>\n");
  exit(1);
}


void mctemp_usage() {
 fprintf(stderr, "--mctemp option usage: --mctemp <temperature in KT>\n");
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
 fprintf(stderr, "--gravity option usage: --gravity <1,2 or 3>\n");
  exit(1);
}

void rstk_usage() {
 fprintf(stderr, "--rstk option usage: --rstk <gravity constant>\n");
  exit(1);
}

void gridmode_usage() {
 fprintf(stderr, "--gridmode option usage: --gridmode <1 or 2>\n");
  exit(1);
}

void vmax_usage() {
 fprintf(stderr, "--vmax option usage: --vmax <maximum number of steps>\n");
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

void em_usage() {
 fprintf(stderr, "--em option usage: --em <EM map definition file>\n");
  exit(1);
}

void proxlim_usage() {
 fprintf(stderr, "--proxlim option usage: --proxlim <distance squared>\n");
  exit(1);
}

void proxmax_usage() {
 fprintf(stderr, "--proxmax option usage: --proxmax <distance squared>\n");
  exit(1);
}

void proxmaxtype_usage() {
 fprintf(stderr, "--proxmaxtype option usage: --proxmaxtype <number of types>\n");
  exit(1);
}


void parse_options(int ministatehandle, int cartstatehandle, int nlig, int argc, char *argv[]) {
  MiniState &ms = ministate_get(ministatehandle);
  CartState &c = cartstate_get(cartstatehandle);
  bool gridspecify = 0;
  for (int n = 0; n < argc; n++) {
    char *arg = argv[n];
    if (!strcmp(arg,"--mc")) {
      ms.imc = 1;
    }
    else if (!strcmp(arg,"--mctemp")) {
      if (argc-n < 2) mctemp_usage();    
      double mctemp = atof(argv[n+1]);
      if (mctemp < 0) mctemp_usage();
      ms.mctemp = mctemp;
      n += 1;
    }    
    else if (!strcmp(arg,"--mcensprob")) {
      if (argc-n < 2) mcensprob_usage();    
      double mcensprob = atof(argv[n+1]);
      if (mcensprob < 0 | mcensprob > 1) mcensprob_usage();
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
    else if (!strncmp(arg,"--traj", 4)) {
      ms.iscore = 2;
    }
    else if (!strcmp(arg,"--fix-receptor")) {
      ms.fixre = 1;
    }
    else if (!strcmp(arg,"--only-rot")) {
      ms.itra = 0;
      ms.iori = 1;
    }
    else if (!strcmp(arg,"--only-trans")) {
      ms.itra = 1;
      ms.iori = 0;
    }    
    else if (!strcmp(arg,"--grid")) {
      gridspecify = 1;
      if (argc-n < 3) grid_usage();
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) grid_usage();
      char *gridf = argv[n+2];
      if (!exists(gridf)) {
        fprintf(stderr, "Grid file %s does not exist\n", gridf);
	grid_usage();
      }
      Grid *g = new Grid;
      g->read(gridf);
      if (g->natoms != c.natom[lig-1]) {
        fprintf(stderr, "Wrong number of atoms for ligand %d:\n  Grid file %s: %d, PDB file: %d\n",lig,gridf,g->natoms,c.natom[lig-1]);
	exit(1);
      }
      g->init_prox(cartstatehandle,ms.proxlim,ms.proxmax,ms.proxmaxtype);
      c.grids[lig-1] = g;
      n += 2;
    }
    else if (!strcmp(arg,"--ens") || (!strcmp(arg,"--ensemble"))) {
      if (argc-n < 3) ens_usage();
      int lig = atoi(argv[n+1]);
      if (lig < 1 || lig > nlig) ens_usage();
      char *ensf = argv[n+2];
      if (!exists(ensf)) {
        fprintf(stderr, "Ensemble file %s does not exist\n", ensf);
	ens_usage();
      }
      read_ens(cartstatehandle,lig-1,ensf,1);
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
      if (gravity <= 0 || gravity > 3) gravity_usage();
      ms.gravity = gravity;
      n += 1;
    }
    else if (!strcmp(arg,"--rstk")) {
      if (argc-n < 2) rstk_usage();    
      double rstk = atof(argv[n+1]);
      ms.rstk = rstk;
      n += 1;
    }    
    else if (!strcmp(arg,"--gridmode")) {
      if (argc-n < 2) gridmode_usage();    
      int gridmode = atoi(argv[n+1]);
      if (gridmode < 1 || gridmode > 2) gridmode_usage();
      ms.gridmode = gridmode;
      n += 1;
    }
    else if (!strcmp(arg,"--vmax")) {
      if (argc-n < 2) vmax_usage();    
      int vmax = atoi(argv[n+1]);
      if (vmax <= 0) vmax_usage();
      ms.ivmax = vmax;
      n += 1;
    }
    else if (!strcmp(arg,"--proxlim")) {
      if (gridspecify) {
        fprintf(stderr, "proxlim cannot be specified after grid\n");
	exit(1);
      }
      if (argc-n < 2) proxlim_usage();    
      double proxlim = atof(argv[n+1]);
      if (proxlim <= 0) proxlim_usage();
      if (proxlim >= ms.proxmax) {
        fprintf(stderr, "proxlim must be smaller than proxmax\n");
	exit(1);
      }
      ms.proxlim = proxlim;
      n += 1;
    }
    else if (!strcmp(arg,"--proxmax")) {
      if (gridspecify) {
        fprintf(stderr, "proxmax cannot be specified after grid\n");
	exit(1);
      }
      if (argc-n < 2) proxmax_usage();    
      double proxmax = atof(argv[n+1]);
      if (proxmax <= 0) proxmax_usage();
      if (proxmax <= ms.proxlim) {
        fprintf(stderr, "proxmax must be larger than proxlim\n");
	exit(1);
      }
      ms.proxmax = proxmax;
      n += 1;
    }
    else if (!strcmp(arg,"--proxmaxtype")) {
      if (gridspecify) {
        fprintf(stderr, "proxmaxtype cannot be specified after grid\n");
	exit(1);
      }    
      if (argc-n < 2) proxmaxtype_usage();    
      int proxmaxtype = atoi(argv[n+1]);
      if (proxmaxtype <= 0) proxmaxtype_usage();
      ms.proxmaxtype = proxmaxtype;
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
      ms.has_globalenergy = 1; //because of val forces => moderest
    }
    else if (!strcmp(arg,"--em")) {
      if (argc-n < 2) em_usage();    
      char *emf = argv[n+1];
      if (!exists(emf)) {
        fprintf(stderr, "EM definition file %s does not exist\n", emf);
	em_usage();
      }
      read_densitymaps_(emf,strlen(emf));
      ms.has_globalenergy = 1;
      n += 1;
    }
    else {
      fprintf(stderr, "Unknown option %s\n", arg);
      exit(1);
    }
  }
}
