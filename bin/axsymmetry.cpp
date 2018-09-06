//Applies explicit-axis symmetry to a DOF file

//usage: ./axsymmetry structures.dat
// <total number of ligands> <ligand> <symmetry>
// <axis x> <axis y> <axis z>
// <origin x> <origin y> <origin z>


#include "max.h"
#include "axsym.h"
#include <cmath>
#include <cstdio>

extern "C" void print_struc_(
 const int &seed,
 char *label0,
 const double &energy,
 const double *energies,
 const int &nlig,
 const int *ens,
 const double *phi,
 const double *ssi,
 const double *rot,
 const double *xa,
 const double *ya,
 const double *za,
 const coors2 &locrests,
 const double *morph,
 const int *nhm,
 const int *nihm,
 const modes2 &dlig,
 const int *has_locrests,
 int len_label
);

extern "C" FILE *read_dof_init_(const char *f_, int nlig, int &line, double (&pivot)[3][MAXLIG], int &auto_pivot, int &centered_receptor, int &centered_ligands, int f_len);

extern "C" int read_dof_(FILE *fil, int &line, int &nstruc, const char *f_, idof2 &ens, dof2 &phi, dof2 &ssi, dof2 &rot, dof2 &xa, dof2 &ya, dof2 &za, coors2 &locrests, dof2 &morph, modes2 &dlig, const int &nlig, const int *nhm, const int *nihm, const int *nrens0, const int *morphing, const int *has_locrests, int &seed, char *&label, const int &all_labels, int f_len);

static int seed;
static char *label;

#include <cstdio>
#include <cstring>
#include <cstdlib>

void usage() {
  fprintf(stderr, "usage: ./axsymmetry structures.dat\n  <total number of ligands> <ligand> <symmetry>\n  <axis x> <axis y> <axis z>\n  <origin x> <origin y> <origin z>\n");
  exit(1);
}

bool exists(const char *f) {
  FILE *fil = fopen(f, "r");
  if (fil == NULL) return 0;
  else {
    fclose(fil);
    return 1;
  }
}

int main(int argc, char *argv[]) {
  int i;

  int nrens[MAXLIG];
  int nhm[MAXLIG];
  int nihm[MAXLIG];

  int ens[MAXLIG];
  double phi[MAXLIG];
  double ssi[MAXLIG];
  double rot[MAXLIG];
  double xa[MAXLIG];
  double ya[MAXLIG];
  double za[MAXLIG];
  double morph[MAXLIG];
  double dlig[MAXLIG][MAXMODE+MAXINDEXMODE];

  int morphing[MAXLIG];
  memset(morphing, 0, MAXLIG*sizeof(int));

  int nr_symcopies[MAXLIG];
  int symcopies[MAXLIG][24*MAXLIG];
  int nr_symtrans;
  SymTrans symtrans[24*MAXLIG];

  double forcefactor[MAXLIG];
  double forcerotation[MAXLIG][9];
  for (int n = 0; n < MAXLIG; n++) {
    double *fr = forcerotation[n];
    fr[0] = 1;
    fr[4] = 1;
    fr[8] = 1;
  }


  for (int n = 0; n < MAXLIG; n++) {
    nhm[n] = 0;
    nihm[n] = 0;
    nrens[n] = 0;
  }

  coors2 locrests;
  int has_locrests[MAXLIG];
  memset(has_locrests, 0, MAXLIG*sizeof(int));

  while (1) {
    if (strcmp(argv[argc-3],"--ens")) break;
    int lignr = atoi(argv[argc-2]) - 1;
    nrens[lignr] = atoi(argv[argc-1]);
    argc -= 3;
  }

  if ((argc-3) % 8) {
    fprintf(stderr, "Wrong number of arguments\n"); usage();
  }
  if (!exists(argv[1])) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    exit(1);
  }

  int nlig0 = atoi(argv[2]);
  for (int n = 0; n < nlig0; n++) {
    nr_symcopies[n] = 1;
    symcopies[n][0] = n;
  }

  AxSymmetry syms[MAXLIG];
  int nr_syms = (argc-3)/8;
  for (int n = 0; n < nr_syms; n++) {
    AxSymmetry &sym = syms[n];
    sym.ligand = atoi(argv[8*n+3]);
    sym.symtype = atoi(argv[8*n+4]);
    sym.angle = 0;
    sym.axis[0] = atof(argv[8*n+5]);
    sym.axis[1] = atof(argv[8*n+6]);
    sym.axis[2] = atof(argv[8*n+7]);
    sym.origin[0] = atof(argv[8*n+8]);
    sym.origin[1] = atof(argv[8*n+9]);
    sym.origin[2] = atof(argv[8*n+10]);
  }

  //read DOFs and set pivots
  double pivot[3][MAXLIG];
  memset(pivot,0,sizeof(pivot));

  int auto_pivot, centered_receptor, centered_ligands;
  int line;
  FILE *fil = read_dof_init_(argv[1], nlig0, line, pivot, auto_pivot, centered_receptor, centered_ligands, strlen(argv[1]));
  if (centered_receptor != centered_ligands) {
    fprintf(stderr, "Receptor and ligands must be both centered or both uncentered\n");
    exit(1);
  }
  if (!(centered_ligands) && auto_pivot) {
    fprintf(stderr, "With uncentered ligands, pivots must be supplied\n");
    exit(1);
  }


  //main loop
  int nstruc = 0;
  int nlig = prepare_axsym_dof(
   nr_syms,
   syms,
   nlig0,
   nhm,
   nihm,
   nrens,
   morphing,
   has_locrests,

   nr_symcopies,
   symcopies,
   nr_symtrans,
   symtrans,

   forcefactor,
   forcerotation
  );

  /*
  for (int i = 0; i < nlig; i++) {
    printf("%d: forcefactor: %.4f\n", i, forcefactor[i]);
    double *m = forcerotation[i];
    printf("%.4f %.4f %.4f\n", m[0], m[1], m[2]);
    printf("%.4f %.4f %.4f\n", m[3], m[4], m[5]);
    printf("%.4f %.4f %.4f\n\n", m[6], m[7], m[8]);
  }
  */
  double newpivot[3][MAXLIG];
  for (int i = 0; i < nlig0; i++) {
    for (int ii = 0; ii < nr_symcopies[i]; ii++)  {
      int target = symcopies[i][ii];
      newpivot[0][target] = pivot[0][i];
      newpivot[1][target] = pivot[1][i];
      newpivot[2][target] = pivot[2][i];
    }
  }

  if (auto_pivot) printf("#pivot auto\n");
  else {
    for (i = 0; i < nlig; i++) {
      printf("#pivot %d %.3f %.3f %.3f\n",
	i+1, newpivot[0][i], newpivot[1][i], newpivot[2][i]);
    }
  }
  if (centered_receptor) printf("#centered receptor: true\n");
  else printf("#centered receptor: false\n");
  if (centered_ligands) printf("#centered ligands: true\n");
  else printf("#centered ligands: false\n");

  while (1) {

    int result = read_dof_(fil, line, nstruc, argv[1], ens, phi, ssi, rot,
     xa, ya, za, locrests,
     morph, dlig, nlig0, nhm, nihm, nrens, morphing, has_locrests,
     seed, label, 1, strlen(argv[1])
    );
    if (result != 0) break;
    apply_axsym(
     nr_symtrans,
     symtrans,
     morph,
     ens,
     phi,
     ssi,
     rot,
     xa,
     ya,
     za,
     dlig,
     has_locrests,
     locrests
    );

    double dummy = 0;
    int lablen = 0;
    if (label) lablen = strlen(label);

    print_struc_(
     seed,
     label,
     dummy,
     NULL,
     nlig,
     ens,
     phi,
     ssi,
     rot,
     xa,
     ya,
     za,
     locrests,
     morph,
     nhm,
     nihm,
     dlig,
     has_locrests,
     lablen
    );
  }
}
