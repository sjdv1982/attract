#combines all ligands in a .dat file with all ligands
import sys
from math import *
from _read_struc import read_struc
import itertools

header,structures = read_struc(sys.argv[1])
structures = list(structures)
ligands = len(structures[0][1])
combinations = long(len(structures))**ligands
if combinations > 10**9:
  raise ValueError("crosscombine.py: more than 1 billion combinations, I give up...")

combi = [list(range(len(structures)))] * ligands

  
stnr = 0
for h in header: print h
for c in itertools.product(*combi): 
  stnr += 1
  print "#"+str(stnr)
  print "## crosscombine",
  for v in c: print v+1,
  print
  s0 = structures[c[0]]
  for l in s0[0]: print l
  for vnr,v in enumerate(c):
    print structures[v][1][vnr]
  
