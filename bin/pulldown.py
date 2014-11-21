"""
Prints out all structure indices in a DAT file that are similar to the ligand coordinates in a provided file
usage: python pulldown.py  <DAT file> <.npy coordinate file> <PDB file> <ligand RMSD cutoff> [collect options]
"""
import sys

import numpy
import collectlibpy as collectlib
from scipy.spatial.distance import cdist
  
import os

dat = sys.argv[1]
refecoorfile = sys.argv[2]
pdb = sys.argv[3]
cutoff = float(sys.argv[4])

initargs = [dat, pdb, pdb] + sys.argv[5:]
collectlib.collect_init(initargs)

#load reference coordinates
refecoor = numpy.load(refecoorfile)
assert len(refecoor.shape) == 2, refecoor.shape
coorsize = refecoor.shape[1]
assert coorsize % 3 == 0, coorsize
natom = coorsize / 3
lim = cutoff * cutoff * natom

strucnr = 0
while 1:
  strucnr += 1
  if not strucnr % 1000: print >> sys.stderr, strucnr
  result = collectlib.collect_next()
  if result: break
  coor = collectlib.collect_coor_raw()    
  if strucnr == 1:
    assert len(coor) == 2 * coorsize, (len(coor), coorsize)
  
  coor = coor[coorsize:]  
  coor2 = coor[numpy.newaxis]
  d = cdist(coor2, refecoor, 'sqeuclidean')[0]   
  if d.min() <= lim:
    print strucnr