import sys

import numpy as np
from math import *
import rmsdlib
from multiprocessing import Pool, Queue

def apply_matrix(atoms, pivot, rotmat, trans):
  ret = []
  for atom in atoms:
    a = atom-pivot
    atom2 = a.dot(rotmat) + pivot + trans
    ret.append(atom2)
  return ret

def write_pdb(outputfile, lines, atoms, extralines):
  outp = open(outputfile, "w")
  count = 0
  pos = 0
  data = zip(lines, atoms)
  while 1:
    while pos < len(extralines):
      p,d = extralines[pos]
      if count < p: break
      print >> outp, d.rstrip("\n")
      pos += 1
    if count == len(data): break
    l,a = data[count]
    ll = l[:30] + "%8.3f%8.3f%8.3f" % (a[0],a[1],a[2]) + l[54:].rstrip("\n")
    print >> outp, ll
    count += 1
  outp.close()

import sys
import argparse
a = argparse.ArgumentParser(prog="rmsd-matrix-pdb.py")
a.add_argument("pdblist")
a.add_argument("--np",type=int)
a.add_argument("--allatoms", action="store_true")
a.add_argument("--ca", action="store_true")
args = a.parse_args()

#structures
pdbfiles = [l.strip().strip("\n") for l in open(args.pdblist) if len(l.strip().strip("\n"))]

#read atoms
pdbs = [rmsdlib.read_pdb(f) for f in pdbfiles]
for pdb in pdbs:
    rmsdlib.check_pdbs((pdb,), pdbs)
coors = np.array([list(pdb.coordinates()) for pdb in pdbs])

#select backbone
if args.allatoms:
  assert not args.ca
else:
  if args.ca:
    atomnames = ("CA",)
  else:
    atomnames = ("CA","C","O","N")
  amask = np.array([a.name in atomnames for a in pdbs[0].atoms()])
  coors = coors[:,amask]

def run(index):
  a2 = coors[index]
  array_a1 = coors[index+1:]
  rotmat, offset, rmsd = rmsdlib.multifit(array_a1, a2)
  return rmsd

pool = Pool(args.np)
try:
  result = pool.map(run, range(len(coors)-1))
except KeyboardInterrupt:
  pool.terminate()
  sys.exit(1)

for n1 in range(len(coors)-1):
  offset = n1 + 1
  for index in range(len(result[n1])):
    n2 = offset + index
    print "%d %d %.3f" % (n1+1, n2+1, result[n1][index])
