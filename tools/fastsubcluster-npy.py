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

npyfile=sys.argv[1]
clustfile=sys.argv[2]
cutoff=float(sys.argv[3])
output_subclust=sys.argv[4]
output_superclust=sys.argv[5]

import numpy as np


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

coors = np.load(sys.argv[1])
assert coors.ndim in (2,3)
if coors.ndim == 2:
    assert not coors.shape[1] % 3, coors.shape
    coors = coors.reshape(len(coors), coors.shape[1]//3, 3)
nstruc = len(coors)
natom = coors.shape[1]
lim = cutoff * cutoff * natom

leafclusts = np.zeros(dtype=float,shape=(leafclustmaxsize,natom, 3))


for rootclustnr, rootclust in enumerate(rootclusters):
    assert min(rootclust) > 0 and max(rootclust) <= nstruc, "Wrong cluster %d" % (rootclustnr+1)

nsubclust=0
for rootclustnr, rootclust in enumerate(rootclusters):
    nleafclust=0
    print >> sys.stderr, rootclustnr+1
    if len(rootclust) == 1:
        subclust.append(rootclust)
        nsubclust += 1
        superclust.append([nsubclust])
        continue
    superclust.append([])
    leafclust = [] #subclusters in the CURRENT rootcluster
    for struc in rootclust:

        coor = coors[struc-1]

        if nleafclust > 0:
            dif = leafclusts[:nleafclust] - coor
            d = np.einsum("ijk,ijk->i", dif, dif)
            if d.min() < lim:
                clus = d.argmin()
                leafclust[clus].append(struc)
                continue

        #create new subcluster
        leafclusts[nleafclust,:] = coor
        leafclust.append([struc])
        nleafclust += 1
        nsubclust+=1
        superclust[rootclustnr].append(nsubclust)
    subclust += leafclust

write_clustfile(superclust, output_superclust)
write_clustfile(subclust, output_subclust)
