import sys
from math import *
from _read_struc import read_struc
try:
  import psyco
  psyco.full()
except:
  pass  

header,structures = read_struc(sys.argv[1])

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

dofs = []
for l1,l2 in structures:
  d = [[],[]]
  for ll in l2:
    d0 = [float(v) for v in ll.split()]    
    d[0].append(euler2rotmat(*d0[:3]))
    d[1] += d0[3:]
  dofs.append(d)

for h in header: print h

radgyr = 30.0 #radius of gyration
ncomp = len(dofs[0][0]) - 1 #ligand RMSD = ignore receptor

stnr = 0
st2nr = 0
lim = 0.02
okdofs = []
for d,s in zip(dofs, structures):
  stnr += 1
  ok = True
  for n,dd in okdofs:
    rmsd = 0
    for rm1, rm2 in zip(d[0],dd[0])[1:]:  #assume fixed receptor...
      for i in range(3):
        r1,r2 = rm1[i],rm2[i]
	dr = (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2])
	dissq = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]
        rmsd += radgyr/ncomp * dissq/3
	if rmsd > lim: break
      if rmsd > lim: break
    if rmsd < lim:
      for p1,p2 in zip(d[1],dd[1]):
        rmsd += (p1-p2)*(p1-p2)/3
        if rmsd > lim: break
    if rmsd < lim: 
      ok = False
      break
  if ok == False: continue
  okdofs.insert(0,(stnr,d))
  st2nr += 1
  print "#"+str(st2nr)
  print "##"+str(stnr) + " => deredundant"
  l1,l2 = s
  for l in l1: print l
  for l in l2: print l
  
