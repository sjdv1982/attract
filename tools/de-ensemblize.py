import sys
from math import *
from _read_struc import read_struc
import random

if len(sys.argv) != 3:
  print >> sys.stderr, "Usage: de-ensemblize <DOF file> <molecule ID> "
  sys.exit()

header,structures = read_struc(sys.argv[1])
molecule = int(sys.argv[2])
assert molecule > 0

stnr = 0
for h in header: print h
for s in structures: 
  stnr += 1
  print "#"+str(stnr)
  l1,l2 = s
  assert len(l2) >= molecule, "Number of molecules is less than %d" % molecule
  head = []
  for l in l2: head.append("")
  head[molecule-1] = str(e) + " "
  for l in l1: print l
  for lnr, l in enumerate(l2): 
    if lnr+1 == molecule: 
      ll = l.split()
      print "    ",
      for f in ll[1:]: print f,
      print
    else:
      print l
