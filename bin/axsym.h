#include "state.h"

extern void axsym_fold_grads(
 MiniState &m, CartState &c,
 double *grads0, double *grads, double *morphgrads
);
extern void prepare_axsym_cartstate(CartState &c);

extern int prepare_axsym_dof(
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
);

extern void apply_axsym(
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
);

extern "C" void euler2rotmat_(const double &phi,const double &ssi, const double &rot, double (&rotmat)[9]);
extern "C" void matinv_(double *rotmat);
extern "C" void matmult_(double *rotmat1, double *rotmat2, double *rotmat);
extern "C" void vecmatmult_(double *v0, double *m, double *v);

