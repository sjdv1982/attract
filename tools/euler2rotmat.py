from math import *

def euler2rotmat(phi,ssi,rot):
  cs=cos(ssi)
  cp=cos(phi)
  ss=sin(ssi)
  sp=sin(phi)
  cscp=cs*cp
  cssp=cs*sp
  sscp=ss*cp
  sssp=ss*sp
  crot=cos(rot)
  srot=sin(rot)

  r1 = crot * cscp + srot * sp  
  r2 = srot * cscp - crot * sp
  r3 = sscp

  r4 = crot * cssp - srot * cp
  r5 = srot * cssp + crot * cp
  r6 = sssp

  r7 = -crot * ss
  r8 = -srot * ss
  r9 = cs
  return ((r1,r2,r3),(r4,r5,r6),(r7,r8,r9))
