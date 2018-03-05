typedef struct {
  int ligand;
  int targetligand;
  double rotmatsym[9];
  double origin[3];
} SymTrans;

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
);
