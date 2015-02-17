import sys, os
import rmsdlib

pdblist = sys.argv[1]
clustfile=sys.argv[2]
cutoff=float(sys.argv[3])
output_subclust=sys.argv[4]
output_superclust=sys.argv[5]

pdbfiles = [l.strip().strip("\n") for l in open(pdblist)]
for pdb in pdbfiles:
  assert os.path.exists(pdb), pdb

import numpy
from scipy.spatial.distance import pdist, squareform
from _read_struc import read_struc
from math import sqrt

def read_clustfile(clustfile):
  clust = []
  for l in open(clustfile):
    ll = l.split()[3:]
    clust.append([int(v) for v in ll])
  return clust  


def write_clustfile(clust, clustfile):
  cf = open(clustfile, "w")
  for cnr, c in enumerate(clust):
    print >> cf, "Cluster %d ->" % (cnr+1), 
    for cc in c: print >> cf, cc,
    print >> cf, ""

rootclusters = read_clustfile(clustfile)

superclust = []
subclust = []

maxstruc = 100000

coor = rmsdlib.read_pdb(pdbfiles[0]).coordinates()
coor = numpy.array(coor)
natom = len(coor)
lim = cutoff * cutoff * natom
coor = coor.flatten()
coorsize = len(coor)
clust_struc = numpy.zeros(dtype=float,shape=(maxstruc,coorsize))    

for rootclustnr, rootclust in enumerate(rootclusters):
  print >> sys.stderr, rootclustnr+1
  if len(rootclust) == 1:
    subclust.append(rootclust)
    superclust.append([len(subclust)])
    continue
  leafclust = []
  for cnr, c in enumerate(rootclust):
    pdb = pdbfiles[c-1]
    coor = rmsdlib.read_pdb(pdb).coordinates()
    coor = numpy.array(coor).flatten()
    clust_struc[cnr,:] = coor
  d = squareform(pdist(clust_struc[:len(rootclust)], 'sqeuclidean'))
  d2 = d<lim
  
  clustered = 0
  while clustered < len(rootclust):
    neigh = d2.sum(axis=0)
    heart = neigh.argmax()
    leaf = numpy.where(d2[heart])[0]    
    for cs in leaf:
      d2[cs,:] = False
      d2[:, cs] = False
    leaf = [heart+1] + [v+1 for v in leaf if v != heart]  
    leafclust.append(leaf)
    clustered += len(leaf)

  mapped_root = []
  for leaf in leafclust:
    mapped_leaf = [rootclust[n-1] for n in leaf]
    subclust.append(mapped_leaf)
    mapped_root.append(len(subclust))
  superclust.append(mapped_root)  

write_clustfile(superclust, output_superclust)
write_clustfile(subclust, output_subclust)
