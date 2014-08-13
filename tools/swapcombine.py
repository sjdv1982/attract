#combines all ligands in a .dat file with all ligands, by taking N-1 ligands from one structure and 1 ligand from another
import sys
from math import *
from _read_struc import read_struc
import itertools

header,structures = read_struc(sys.argv[1])
structures = list(structures)
ligands = len(structures[0][1])
assert ligands > 1
combinations = long(len(structures))**2 * ligands
if combinations > 10**9:
  raise ValueError("swapcombine.py: more than 1 billion combinations, I give up...")

combi = [list(range(len(structures)))] * ligands

  
stnr = 0
for h in header: print h
for varlig in reversed(range(ligands)):
  for n in range(len(structures)):
    x0, x1 = structures[n]
    for nn in range(len(structures)):
      if n == nn and varlig != ligands-1: continue
      y0, y1 = structures[nn]
      stnr += 1
      print "#"+str(stnr)
      print "## swapcombine",
      for lnr in range(ligands): 
        v = nn+1 if lnr==varlig else n +1
        print v,
      print  
      for l in x0: print l
      for lnr in range(ligands):
        print y1[lnr] if lnr==varlig else x1[lnr]
      
