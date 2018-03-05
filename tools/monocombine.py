#replaces ligand N in .dat file 1 with ligand N in .dat file 2. .dat file 2 must contain a single structure
import sys
from math import *
from _read_struc import read_struc
import itertools

ligand = int(sys.argv[1])
header,structures2 = read_struc(sys.argv[3])
structures2 = list(structures2)
assert len(structures2) == 1, sys.argv[3]
other = structures2[0]
assert ligand <= len(other[1])

header,structures = read_struc(sys.argv[2])
  
stnr = 0
for h in header: print h
for s in structures: 
  assert len(s[1]) == len(other[1])
  stnr += 1
  print "#"+str(stnr)
  for l in s[0]: print l
  for vnr,v in enumerate(s[1]):
    if vnr +1 == ligand:
      print other[1][vnr]
    else:
      print s[1][vnr]
  
