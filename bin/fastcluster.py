"""
Fast clustering of two-body docking solutions
usage: python fastcluster.py <DAT file> <PDB file> <cutoff> [collect options]

"""
import sys

import numpy
import collectlibpy as collectlib
from scipy.spatial.distance import cdist
  
import os

dat = sys.argv[1]
pdb = sys.argv[2]
cutoff = float(sys.argv[3])

initargs = [dat, pdb, pdb] + sys.argv[4:]
collectlib.collect_init(initargs)

first = True
clustmaxsize = 1000

strucnr = 1
result = collectlib.collect_next()
coor = collectlib.collect_coor_raw()    
assert len(coor) % 6 == 0, len(coor)
coorsize = len(coor) / 2
  
natom = len(coor) / 6
lim = cutoff * cutoff * natom

clusts = numpy.zeros(dtype=float,shape=(clustmaxsize,coorsize))    
coor = coor[coorsize:]
clusts[0,:] = coor
nclust = 1
first = False        
clustind = [[1]]

while 1:
  strucnr += 1
  result = collectlib.collect_next()
  if result: break
  coor = collectlib.collect_coor_raw()    
  if nclust == clustmaxsize:
    clustmaxsize = int(clustmaxsize*1.2)
    clusts2 = numpy.zeros(dtype=float,shape=(clustmaxsize,coorsize))    
    clusts2[:nclust] = clusts[:]
    clusts = clusts2      
  
  coor = coor[coorsize:]  
  coor2 = coor[numpy.newaxis]

  d = cdist(coor2, clusts[:nclust], 'sqeuclidean')[0]   
  if d.min() < lim:
    c = d.argmin()
    clustind[c].append(strucnr)
    continue
    
  #create new cluster  
  clusts[nclust,:] = coor
  clustind.append([strucnr])
  nclust += 1 
  print >> sys.stderr, strucnr, nclust

for cnr, c in enumerate(clustind):
  print "Cluster", cnr+1, "->", 
  for cc in c: print cc,
  print