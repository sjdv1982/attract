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
    receptorid = []
    data = open(pdb).readlines()
    data = [x for x in data if 'ATOM' in x]
    for count,line in enumerate(data):
        tmp = line.replace('-',' -')
        list = tmp.split()
        receptor.append((count+1,int(list[4]),float(list[5]),float(list[6]),float(list[7])))
                          
    return receptor 

  
def make(pdb,nlistcut=30):        
    atomlist = read_pdb(pdb)
    coor = [[x,y,z] for id, res, x, y, z in atomlist]
    atomlist = [atom[0] for atom in atomlist]
    coor = np.array(coor)
    Y = scipy.spatial.cKDTree(coor)
    nlist = []
    #search for neighbors and bonds
    for atom in atomlist:
      dist, neighbor = Y.query(coor[atom-1:atom],k=nlistcut+1)
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