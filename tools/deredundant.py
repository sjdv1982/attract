import sys
from math import *
from _read_struc import read_struc
try:
  import psyco
  psyco.full()
except:
  pass  

header,structures = read_struc(sys.argv[1])
structures = list(structures)
from euler2rotmat import euler2rotmat

dofs = []
for l1,l2 in structures:
  d = [[],[]]
  for ll in l2:
    d0 = [float(v) for v in ll.split()]    
    d[0].append(euler2rotmat(*d0[:3]))
    d[1].append(d0[3:])
  dofs.append(d)

for h in header: print h

radgyr = 30.0 #radius of gyration
ncomp = len(dofs[0][0]) - 1 #ligand RMSD = ignore receptor

stnr = 0
st2nr = 0
lim = 0.05
if len(sys.argv) > 2:
   lim = float(sys.argv[2])

okdofs = []
for d,s in zip(dofs, structures):
  stnr += 1
  ok = True
  okdofnr = 0
  for n,dd in okdofs:
    okdofnr += 1
    rmsd = 0
    for lignr in range(1, len(d[0])):  #assume fixed receptor...
      rm1 = d[0][lignr], d[1][lignr]
      rm2 = dd[0][lignr], dd[1][lignr]
      for i in range(3):
        r1,r2 = rm1[0][i],rm2[0][i]        
	dr = (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2])
	dissq = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]
        rmsd += radgyr/ncomp * dissq/3
	if rmsd > lim: break
      if rmsd > lim: break
      
      r1,r2 = rm1[1],rm2[1]
      dr = (r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2])
      dissq = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]
      rmsd += dissq/ncomp
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
  
