import sys, os, tempfile

receptor = "/dev/null"
anr = 0
while 1:
  anr += 1
      
  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]

  if anr <= len(sys.argv)-2 and arg == "--receptor":
    receptor = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

datfile=sys.argv[1]
clustfile=sys.argv[2]
ligandpdbfile=sys.argv[3]
cutoff=float(sys.argv[4])
output_subclust=sys.argv[5]
output_superclust=sys.argv[6]
rest = sys.argv[7:]

import numpy
import collectlibpy as collectlib
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

tmpfd, tmpnam = tempfile.mkstemp(text=True)
tmpf = os.fdopen(tmpfd, "w")
try:
  header, structures = read_struc(datfile)
  structures = list(structures)
  for h in header: print >> tmpf, h
  snr = 0
  for c in rootclusters:
    for cc in c:
      snr += 1
      print >> tmpf, "#%d" % snr
      s1, s2 = structures[cc-1]
      for l in s2: print >> tmpf, l      
  tmpf.close()

  assert snr == sum([len(c) for c in rootclusters]), (snr, sum([len(c) for c in rootclusters]))

  superclust = []
  subclust = []

  initargs = [tmpnam, receptor, ligandpdbfile] + rest
  maxstruc = 100000

  collectlib.collect_init(initargs)
  result = collectlib.collect_next()
  coor = collectlib.collect_coor_raw()
  coorstart = collectlib.ieins[0]
  assert len(coor) % 3 == 0, len(coor)
  coor = coor[3*coorstart:]
  coorsize = len(coor)
  natom = len(coor) / 3
  lim = cutoff * cutoff * natom
  clust_struc = numpy.zeros(dtype=float,shape=(maxstruc,coorsize))    
  counter = 0

  clust_struc[0,:] = coor

  order_err = "Clusters must be ordered: Cluster 1 =>  1 2 3, Cluster 2 => 4 5 6, etc. "
  for rootclustnr, rootclust in enumerate(rootclusters):
    print >> sys.stderr, rootclustnr+1
    if len(rootclust) == 1:      
      assert rootclust[0] == counter+1, order_err + "Expected: %d; read %d" % (counter+1, rootclust[0])
      subclust.append(rootclust)
      superclust.append([len(subclust)])
      if rootclustnr > 0:
        result = collectlib.collect_next()
        counter += 1
      continue
    leafclust = []
    for cnr, c in enumerate(rootclust):
      assert c == counter+1, order_err + "Expected: %d; read %d" % (counter+1, c)
      if cnr == 0 and rootclustnr == 0: 
        counter += 1
        continue
      result = collectlib.collect_next()      
      coor = collectlib.collect_coor_raw()      
      clust_struc[cnr,:] = coor[3*coorstart:]
      counter += 1
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
finally:
  os.remove(tmpnam)
