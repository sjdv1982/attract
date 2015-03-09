import sys, os, tempfile

leafclustmaxsize = 100000

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
from _read_struc import read_struc
from scipy.spatial.distance import cdist

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
  del structures

  assert snr == sum([len(c) for c in rootclusters]), (snr, sum([len(c) for c in rootclusters]))

  superclust = []
  subclust = []

  initargs = [tmpnam, receptor, ligandpdbfile] + rest

  collectlib.collect_init(initargs)
  result = collectlib.collect_next()
  coor = collectlib.collect_coor_raw()
  coorstart = collectlib.ieins[0]
  assert len(coor) % 3 == 0, len(coor)
  coor = coor[3*coorstart:]
  coorsize = len(coor)
  natom = len(coor) / 3
  lim = cutoff * cutoff * natom
  
  leafclusts = numpy.zeros(dtype=float,shape=(leafclustmaxsize,coorsize))    
  
  nsubclust=0
  for rootclustnr, rootclust in enumerate(rootclusters):
    nleafclust=0
    print >> sys.stderr, rootclustnr+1
    if len(rootclust) == 1:
      subclust.append(rootclust)
      superclust.append([nsubclust])
      nsubclust += 1
      if rootclustnr > 0:
        result = collectlib.collect_next()
      continue
    superclust.append([])
    leafclust = [] #subclusters in the CURRENT rootcluster
    for cnr, c in enumerate(rootclust):
      strucnr=rootclust[cnr]
      if cnr == 0 and rootclustnr == 0:
        pass
      else:   
        result = collectlib.collect_next()
        coor = collectlib.collect_coor_raw()
        coor = coor[3*coorstart:]
      
      coor2 = coor[numpy.newaxis]
      
      if nleafclust > 0:
        d = cdist(coor2, leafclusts[:nleafclust], 'sqeuclidean')[0]   
        if d.min() < lim:
          c = d.argmin()
          leafclust[c].append(strucnr)
          continue
    
      #create new subcluster  
      leafclusts[nleafclust,:] = coor
      leafclust.append([strucnr])
      nleafclust += 1
      nsubclust+=1
      superclust[rootclustnr].append(nsubclust)      
    subclust += leafclust  

  write_clustfile(superclust, output_superclust)
  write_clustfile(subclust, output_subclust)
finally:
  os.remove(tmpnam)
