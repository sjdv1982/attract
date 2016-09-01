import sys,os
import numpy as np

name = sys.argv[1]
cutoff = sys.argv[2]

os.system('rm final.pdb')
os.system('$ATTRACTTOOLS/top '+name+'.dat 5000 > '+name+'-top5000.dat')
os.system('$ATTRACTDIR/fix_receptor '+name+'-top5000.dat 2 > '+name+'-top5000.dat-fixre')
os.system('$ATTRACTTOOLS/backbone ubB-aa.pdb > ubB-bb.pdb')
os.system('$ATTRACTDIR/matrix-lrmsd '+name+'-top5000.dat-fixre ubA-aa.pdb ubB-bb.pdb > '+name+'-top5000.matrix-lrmsd')

os.system('$ATTRACTDIR/cluster_struc '+name+'-top5000.matrix-lrmsd '+cutoff+' 1 > '+name+'-top5000.clust'+cutoff)

ene = open(name+'.dat').readlines()
ene = [float(l.split()[-1]) for l in ene if 'Energy' in l]
cluster = []
for l in open(name+'-top5000.clust'+cutoff).readlines():
  tmp = l.split()[3:]
  tmp = [int(t) for t in tmp]
  tmp.sort()
  totalene = 0.0#ene[tmp[0]-1]
  beststruc = tmp[0]
  for t in tmp[:4]:
    totalene += ene[t-1]
      
  totalene /= float(len(tmp[:4]))
  cluster.append((totalene,beststruc))
  
cluster.sort()

out = open(name+'-top5000-clust'+cutoff+'.score','w')
for i, item in enumerate(cluster):
  out.write(str(i+1)+' '+str(item[0])+' '+str(item[1])+'\n')
  
out.close()
os.system('python $ATTRACTTOOLS/select-structures.py '+name+'-top5000.dat-fixre '+" ".join([str(item[1]) for item in cluster[:100]])+' > '+name+'-top5000-top100cluster.dat')
os.system('$ATTRACTDIR/collect '+name+'-top5000-top100cluster.dat ubA-aa.pdb ubB-aa.pdb > '+name+'-top5000-top100cluster.pdb')
os.system('ln -s '+name+'-top5000-top100cluster.pdb final.pdb')