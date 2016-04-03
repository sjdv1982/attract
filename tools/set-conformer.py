import sys
from math import *
from _read_struc import read_struc
import random

if len(sys.argv) != 4:
  print >> sys.stderr, "Usage: set-conformer <DOF file> <molecule ID> <conformer>"
  sys.exit()

header,structures = read_struc(sys.argv[1])
molecule = int(sys.argv[2])
assert molecule > 0
conformer = int(sys.argv[3])
assert conformer > 0

stnr = 0
for h in header: print h
for s in structures: 
  stnr += 1
  print "#"+str(stnr)
  l1,l2 = s
  assert len(l2) >= molecule, "Number of molecules is less than %d" % molecule
  for l in l1: print l
  for lnr,l in enumerate(l2):
    if lnr == molecule-1:
      ll = l.split()
      if len(ll) == 7:
        print str(conformer), " ".join(ll[1:])
      else:
        print str(conformer), l
    else:
      print l
  
