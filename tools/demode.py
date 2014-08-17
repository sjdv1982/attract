"""
Usage: demode <structure file> [ligand]
Removes all modes from a ligand
if ligand is not specified, then remove the modes from all ligands
"""
import sys
from _read_struc import read_struc
import random
from math import *
header,structures = read_struc(sys.argv[1])
ligand = None
ens = [0,0]
if '--ens' in sys.argv:
  i = sys.argv.index('--ens')
  ensligand = int(sys.argv[i+1])
  ens[ensligand-1] = 1
  sys.argv = sys.argv[:i]+sys.argv[i+2:]
  
if len(sys.argv) > 2:
  ligand = int(sys.argv[2])
  
for h in header: print h
stnr = 0
for s in structures:
  stnr += 1
  l1,l2 = s
  for lnr in range(len(l2)):
    if ligand is not None and lnr+1 != ligand: continue
    l = l2[lnr]
    values = [float(v) for v in l.split()]  
    if len(values) > 6+ens[lnr]:
      values = values[:6+ens[lnr]]
      l2[lnr] = "  " + " ".join([("%.6f" % v) for v in values]) 
  print "#"+str(stnr)
  for l in l1: print l
  for l in l2: print l

