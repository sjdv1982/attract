"""
  program  modes
 
  this program calculates the Hinsen modes (Proteins, 1998) using only CA atoms of a protein
  the modes are adapted to the full protein structure by moving each residue as rigid body
  usage: $path/modes structure.pdb
  author: Martin Zacharias, Jacosb University Bremen
  
  Ported to Python by Sjoerd de Vries, 2014
"""

import numpy, sys
from math import *

def read_pdb(f):
  coor, res, atom = [], [], []
  curr_resid = None
  resindex = 0
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    coor.append((x,y,z))
    resid = l[21:26]
    if resid != curr_resid: 
      curr_resid = resid
      resindex += 1      
    res.append(resindex)  
    atom.append(l[12:16])
  return coor, res, atom

  
coor, res, atom = read_pdb(sys.argv[1])  

ca = []
for ccoor, catom in zip(coor, atom):
  if catom != " CA ": continue
  ca.append(numpy.array(ccoor))
asize = 3*len(ca)  
d2f = numpy.zeros((asize,asize))
for i in range(len(ca)):
  for j in range(len(ca)):
    if i == j: continue
    dx = ca[i] - ca[j]
    r2 = dx.dot(dx)
    fb = exp(-r2/16.0)
    b0 = 2 * fb/r2
    ii = 3 * i
    jj = 3 * j
    for k in range(3):
      for l in range(3):
        d2f[ii+k,jj+l] -= b0*dx[k]*dx[l] 
    for k in range(3):
      for l in range(k,3):
        d2f[ii+k,ii+l] += b0*dx[k]*dx[l] 

for i in range(asize):
  for j in range(0, i+1): 
    d2f[i,j] = d2f[j,i] 

U,s,V = numpy.linalg.svd(d2f)    
ind = numpy.argsort(s)[::-1]
U = U[:, ind]
s = s[ind]
V = V[:, ind]
for i in range(20):  
  n = asize-7-i
  vec = [-float(x) for x in V[n]]
  print i+1, float(s[n])
  for resindex in res:
    pos = 3 * (resindex - 1)
    for n in range(3):
      print vec[pos+n],
      
  print
