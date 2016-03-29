#!/usr/bin/env python

import subprocess
import numpy as np
import sys, os

name = 'result-amber-sorted'
  
benchmark = False
if '--benchmark' in sys.argv:
  benchmark = True

subprocess.call('$ATTRACTTOOLS/splitmodel '+name+'.pdb > '+name+'.modellist',shell=True)
subprocess.call('$FCCHOME/make_contacts.py -f '+name+'.modellist -n 4',shell=True)
subprocess.call("sed -e 's/pdb/contacts/' "+name+".modellist | sed -e '/^$/d' > "+name+".contactlist",shell=True)
subprocess.call('$FCCHOME/calc_fcc_matrix.py -f '+name+'.contactlist -o '+name+'.fcc_matrix',shell=True)
subprocess.call('$FCCHOME/cluster_fcc.py '+name+'.fcc_matrix 0.6 -o '+name+'.fcc_cluster',shell=True)

if benchmark:
  irmsd = np.loadtxt(name+'.irmsd',usecols=[1],unpack=True,ndmin=1)
cluster = []
for line in open(name+'.fcc_cluster').readlines():
  tmp = line.split()
  tmpclus = [int(i)-1 for i in tmp[3:]]
  tmpclus.sort()
  cluster.append(tmpclus)
 
energy = np.loadtxt(name+'.score',unpack=True,usecols=[1]) 
keep = []
for j,c in enumerate(cluster):
  if benchmark:
    irmsdc = [(energy[i],i,irmsd[i]) for i in c]
  else:
    irmsdc = [(energy[i],i,100.0) for i in c]

  irmsdc.sort()
  irmsdmin = 100.0
  ene = 0.0
  count = 0
  save = irmsdc[0][1]
  savenr = 1
  for i in range(len(irmsdc)):
    if i == 4:
      break
    ene += irmsdc[i][0]
    count += 1
    if benchmark and irmsdmin > irmsdc[i][2]:
      irmsdmin = irmsdc[i][2]
      save = irmsdc[i][1]
      savenr = i+1
  
  if benchmark:
    keep.append((ene/float(count),save,savenr,irmsdmin))
  else:
    keep.append((ene/float(count),save,savenr))
  
keep.sort()
if benchmark:
  os.system('python $ATTRACTTOOLS/joinmodel.py '+" ".join([name+'-'+str(i+1)+'.pdb' for j,i,k,l in keep])+' > '+name+'-fcc.pdb')
else:  
  os.system('python $ATTRACTTOOLS/joinmodel.py '+" ".join([name+'-'+str(i+1)+'.pdb' for j,i,k in keep])+' > '+name+'-fcc.pdb')
subprocess.call('rm '+name+'-[0-9]*.pdb',shell=True)
out = open(name+'-fcc.score','w')
out.write("# CLUSTERNR CLUSTERENERGY STRUCTURENR CLUSTERMEMBERNR\n") 
for i,item in enumerate(keep):
  out.write(str(i+1)+' '+str(item[0])+' '+str(item[1]+1)+' '+str(item[2])+'\n')
out.close()

if benchmark:
  out2 = open(name+'-fcc.irmsd','w')
  for i,item in enumerate(keep):
    if item[1] <= 2.0:
      print os.path.split(os.getcwd())[1], "cluster", i+1, "near-native, structure", item[1],"cluster member",item[2], item[3]
    out2.write(str(i+1)+' '+str(item[3])+'\n')
  
  out2.close()
