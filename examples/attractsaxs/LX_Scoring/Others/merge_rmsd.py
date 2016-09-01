import numpy as np
import sys

numsym = int(sys.argv[1])
data = sys.argv[2]

datatype = data.split('.')[-1]

if datatype == 'irmsd' or datatype == 'lrmsd':
  rmsds = []
  rmsds.append(np.genfromtxt(data)[:,-1])
  for i in range(numsym):
    rmsds.append(np.genfromtxt(data+'-sym'+str(i+1))[:, -1])
  rmsds = np.array(rmsds).T
  rmsds = np.sort(rmsds, axis = 1)
  rmsds = rmsds.T[0]
  for r,rmsd in enumerate(rmsds):
    print r, rmsd
else:
  rmsds = []
  rmsds.append(np.genfromtxt(data))
  for i in range(numsym):
    rmsds.append(np.genfromtxt(data+'-sym'+str(i+1)))
  rmsds = np.array(rmsds).T
  rmsds = -np.sort(-rmsds, axis = 1)
  rmsds = rmsds.T[0]
  
  for r,rmsd in enumerate(rmsds):
    print rmsd
    
    
