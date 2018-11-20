#!/usr/bin/env python3

"""
fastcluster script
copyright 2016-2017 Sjoerd de Vries, Isaure Chauvot de Beauchene
"""

import numpy as np
import sys, argparse

def npy2to3(npy):
    if len(npy.shape) == 2:
        if npy.shape[1] == 3:
            npy = npy.reshape(1, npy.shape[0], npy.shape[1])
        else:
            npy = npy.reshape(npy.shape[0], npy.shape[1]/3, 3)
    else:
        assert len(npy.shape) == 3
    return npy

def npy3to2(npy):
    if len(npy.shape) == 3:
        npy = npy.reshape(npy.shape[0], 3*npy.shape[1])
    else:
        assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    return npy

def cluster(structures, threshold, already_clustered, chunksize, assign_structures):
    """Clusters structures using an RMSD threshold
    First structure becomes a cluster,
     second structure only if it doesn't cluster with the first, etc.

    structures: 2D numpy array, second dimension = 3 * natoms
      structures must already have been fitted!
    threshold: RMSD threshold (A)
    already_clustered: if nonzero, the first already_clustered structures are
     considered clusters
    chunksize: number of structures to put in a chunk
      This is an implementation detail that only affects the speed, not the result
    """
    if len(structures.shape) == 3:
        assert structures.shape[2] == 3
        structures = structures.reshape(structures.shape[0], structures.shape[1]*3)
    if len(structures.shape) == 2:
        assert structures.shape[1] % 3 == 0
    natoms = structures.shape[1]/3
    # threshold2 = sum-of-sd threshold = (RMSD threshold **2) * natoms
    threshold2 = threshold**2 * natoms

    nclus = 1
    assert already_clustered >= 0
    if already_clustered == 0:
        already_clustered = 1 ## the first structure is always a cluster
    clus_space = 99 + already_clustered
    clus = np.zeros((clus_space, structures.shape[1]))
    clus[:already_clustered] = structures[:already_clustered]
    clustids = list(range(already_clustered))
    for n in range(already_clustered, len(structures), chunksize):
        print("{0}/{1}".format(n, len(structures)), file=sys.stderr)
        sys.stderr.flush()
        chunk = structures[n:n+chunksize]
        d = chunk[:, np.newaxis, :] - clus[np.newaxis, :, :]
        inter_sd = np.einsum("...ij,...ij->...i", d, d)
        #close_inter is a 2D Boolean matrix:
        #  True  (1): chunk[i] is close to (within RMSD threshold of) clus[j]
        #  False (0): chunk[i] is not close to clus[j]
        close_inter = (inter_sd < threshold2)
        # newclustered contains all structures in the chunk that *don't* cluster with an existing cluster
        newclustered = []
        for chunk_index, closest_inter in enumerate(np.argmax(close_inter,axis=1)):
            # closest_inter contains the *first* index of close_inter
            #   with the highest value of close_inter
            # We are interested in the case where close_inter is all False (=> new cluster)
            # In that case, the highest value of close_inter is False, and closest_inter is 0
            # If close_inter is *not* all False (=> existing cluster), one of these conditions is False
            if closest_inter == 0 and close_inter[chunk_index, 0] == False:
                newclustered.append(chunk_index)

        if len(newclustered):
            # Now we have newclustered: the *chunk* index of all structures in the chunk that will be in new clusters
            # Now we want to cluster them among themselves, and add the *structure* id of each new cluster
            chunk_newclustered = chunk[newclustered]
            d = chunk_newclustered[:, np.newaxis, :] - chunk_newclustered[np.newaxis, :, :]
            intra_sd = np.einsum("...ij,...ij->...i", d, d)
            close_intra = (intra_sd < threshold2)

            # set all upper-triangular indices to False
            close_intra[np.triu_indices(len(chunk_newclustered))] = 0
            for nn in range(len(chunk_newclustered)):
                # same logic as for closest_inter;
                #  except that we don't have the chunk index, but the chunk_newclustered index (nn)
                #  and, since we modify close_intra in the "else" clause, argmax is computed later
                closest_intra = np.argmax(close_intra[nn])
                if closest_intra == 0 and close_intra[nn, 0] == False:
                    chunk_index = newclustered[nn]
                    # if clus is full, re-allocate it as a 50 % larger array
                    if nclus == clus_space:
                        clus_space = int(clus_space*1.5)
                        clus_old = clus
                        clus = np.zeros((clus_space, structures.shape[1]))
                        clus[:nclus] = clus_old
                    clus[nclus] = chunk[chunk_index]
                    clustids.append(n+chunk_index)
                    nclus += 1
                else:  # in addition, if we aren't a new cluster, being close to us doesn't matter
                    close_intra[:, nn] = False

    # After assigning the cluster centers,
    #  assign all structures to the closest cluster
    clusters = {a:[a] for a in clustids}
    if assign_structures:
        for n in range(0, len(structures), chunksize):
            chunk = structures[n:n+chunksize]
            d = chunk[:, np.newaxis, :] - clus[np.newaxis, :, :]
            inter_sd = np.einsum("...ij,...ij->...i", d,d)
            best = np.argmin(inter_sd, axis=1)
            for nn in range(len(chunk)):
                bestclust = clustids[best[nn]]
                if bestclust == (n+nn):
                    continue
                clusters[bestclust].append(n+nn)

    return clusters, clustids

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="np array")
parser.add_argument('cutoff', help="cutoff",type=float)
parser.add_argument("--chunk", type=int)
parser.add_argument("--skip", type=int, help="assume that the first <skip> structures have been clustered already",default=0)
parser.add_argument("--no-assign-structures", action="store_true", help="after initial clustering, do not assign each structure to the closest cluster")
args = parser.parse_args()
############
structures = np.load(args.npy)
print(structures.shape, file=sys.stderr)
threshold = args.cutoff
chunksize = 100
if args.chunk:
    chunksize = args.chunk

clustlist = open(args.npy.split(".npy")[0]+"-clust"+str(threshold), "w")
clustnpy = args.npy.split(".npy")[0]+"-clust"+str(threshold)+".npy"

assign_structures = (not args.no_assign_structures)
c, clustids = cluster(structures, threshold, args.skip, chunksize, assign_structures)
csort = sorted(c, key=lambda k: len(c[k]), reverse=True)
npclust = structures[csort]
for nk, k in enumerate(csort):
    # !!! changed to start at 1 !!! 18/12/17
    print("cluster %i => "%(nk+1), end='', file=clustlist),
    for j in c[k]:
        print("%i "%(j+1), end=' ', file=clustlist)
    print(file=clustlist)

np.save(clustnpy, npclust)
clustlist.close()
