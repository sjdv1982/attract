"""
This script generates the neighbor pairs for the unbound protein structure
This is used to generate the restraints file
"""

import numpy as np
import scipy.spatial
import cPickle as pickle
from multiprocessing import Pool
resultlist = []

  
def read_pdb(pdb):
    # Make list of corresponding atoms
    receptor = []
    data = open(pdb).readlines()
    data = [x for x in data if 'ATOM' in x]
    for count,l in enumerate(data):
      x,y,z = [float(f) for f in (l[30:38],l[38:46],l[46:54])]
      resid = l[21:26]
      receptor.append((count+1,resid,x,y,z))
    return receptor 

  
def make(pdb,nlistcut=30):        
    atomlist = read_pdb(pdb)
    coor = [[x,y,z] for id, res, x, y, z in atomlist]
    atomlist = [atom[0] for atom in atomlist]
    coor = np.array(coor)
    Y = scipy.spatial.cKDTree(coor)
    nlist = []
    #search for neighbors and bonds
    k=min(nlistcut+1, len(atomlist))
    for atom in atomlist:
      dist, neighbor = Y.query(coor[atom-1:atom],k=k)
      dist = dist[0]
      neighbor = neighbor[0]
      rc = float(int(max(dist)+0.99999999))
      dist, neighbor = Y.query(coor[atom-1:atom],k=999, distance_upper_bound=rc)
      dist = dist[0]
      neighbor = neighbor[0]        
      nlist.append([(neighbor[i]+1,dst) for i,dst in enumerate(dist) if not atom == neighbor[i]+1 and not dst == np.inf])

    pickle.dump(nlist,open(pdb+'.tree',"wb"))
    
if __name__ == "__main__":
  import sys
  make(sys.argv[1])