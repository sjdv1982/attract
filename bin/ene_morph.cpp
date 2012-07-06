#include <cstdio>

typedef double Coor[3];

extern "C" void grad_morph_ (
  const Coor *forces,
  const int &natom,
  const double &morph,
  const Coor *cmorphd,
  double &deltamorph
) 
{
  if (morph < 0) return;
  double delta = 0;
  for (int n = 0; n < natom; n++) {
    const Coor &f = forces[n];
    const Coor &d = cmorphd[n];
    delta += f[0]*d[0];
    delta += f[1]*d[1];
    delta += f[2]*d[2];
  }
  deltamorph = delta;
  //printf("X %.3f %.3f\n", morph, deltamorph);
}


extern "C" void ene_morph_ (
  const double &fconstant,
  double *morph,
  double *deltamorph,  
  const int &nlig,
  double &ene
  )
{
  ene = 0;
  for (int n = 0; n < nlig; n++) {
    double cmorph = morph[n];
    if (cmorph < 0) continue;
    deltamorph[n] += -2 * fconstant * cmorph;
    //printf("%d %.3f %.3f\n", n+1, cmorph, deltamorph[n]);    
    ene += fconstant * cmorph * cmorph;
  }

}
