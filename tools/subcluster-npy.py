"""
fastcluster script
copyright 2016-2017 Sjoerd de Vries, Isaure Chauvot de Beauchene
"""

import sys
import numpy as np
from scipy.spatial.distance import pdist, squareform
from math import sqrt

npy_file = sys.argv[1]
clustfile = sys.argv[2]
cutoff = float(sys.argv[3])
output_subclust = npy_file.split(".npy")[0]+"-clust"+str(cutoff)
output_superclust = npy_file.split(".npy")[0]+"-clust"+str(cutoff)+'.npy'

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
maxstruc = 30000

coors = np.load(npy_file)
if len(coors.shape) == 2:    
    coors = coors.reshape(len(coors), coors.shape[1]/3, 3)
nstruc, natom, _ = coors.shape
lim = cutoff * cutoff * natom
coorsize = 3 * natom
coors = coors.reshape(nstruc, coorsize)
clust_struc = np.zeros(dtype=float,shape=(maxstruc,coorsize))

for rootclustnr, rootclust in enumerate(rootclusters):
    print >> sys.stderr, rootclustnr+1
    if len(rootclust) == 1:
        subclust.append(rootclust)
        superclust.append([len(subclust)])
        continue
    elif len(rootclust) > maxstruc:
        raise Exception("Cluster %d has length %d, larger than maxstruc=%d" % (rootclustnr+1,len(rootclust), maxstruc))
    leafclust = []
    for cnr, c in enumerate(rootclust):
        coor = coors[c-1]
        clust_struc[cnr,:] = coor
    d = squareform(pdist(clust_struc[:len(rootclust)], 'sqeuclidean'))
    d2 = d<lim
    del d
    
    clustered = 0
    while clustered < len(rootclust):
        neigh = d2.sum(axis=0)
        heart = neigh.argmax()
        leaf = np.where(d2[heart])[0]
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
