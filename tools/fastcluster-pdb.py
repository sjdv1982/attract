"""
Fast all-atom clustering of single-body PDBs
usage: python fastcluster.py <PDB file list> <cutoff>

"""
import sys
import rmsdlib
import numpy
from scipy.spatial.distance import cdist
  
import os

pdblist = sys.argv[1]
cutoff = float(sys.argv[2])
pdbfiles = [l.strip().strip("\n") for l in open(pdblist)]
for pdb in pdbfiles:
  assert os.path.exists(pdb), pdb

clustmaxsize = 1000

strucnr = 1
coor = rmsdlib.read_pdb(pdbfiles[0]).coordinates()
coor = numpy.array(coor)
natom = len(coor)
lim = cutoff * cutoff * natom
coor = coor.flatten()
coorsize = len(coor)

clusts = numpy.zeros(dtype=float,shape=(clustmaxsize,coorsize))    
clusts[0,:] = coor
nclust = 1
clustind = [[1]]

for pdbnr,pdb in enumerate(pdbfiles[1:]):
  strucnr = pdbnr + 2
  coor = rmsdlib.read_pdb(pdb).coordinates()
  coor = numpy.array(coor).flatten()
  assert len(coor) == coorsize, pdb
  if nclust == clustmaxsize:
    clustmaxsize = int(clustmaxsize*1.2)
    clusts2 = numpy.zeros(dtype=float,shape=(clustmaxsize,coorsize))    
    clusts2[:nclust] = clusts[:]
    clusts = clusts2      
  
  d = cdist(coor[numpy.newaxis], clusts[:nclust], 'sqeuclidean')[0]
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