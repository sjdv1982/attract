"""
Usage: de-ensemblize <structure file> [ligand]
Removes the ensemble conformation from a ligand
if ligand is not specified, then remove the modes from all ligands
"""

import sys
from math import *
from _read_struc import read_struc
import random

assert len(sys.argv) in (2,3)
header,structures = read_struc(sys.argv[1])
ligand = None
if len(sys.argv) == 3:
  ligand = int(sys.argv[2])
  assert ligand > 0

stnr = 0
for h in header: print h
for s in structures: 
  stnr += 1
  print "#"+str(stnr)
  l1,l2 = s
  if ligand is not None:
    assert len(l2) >= ligand, "Number of ligands is less than %d" % ligand
  head = []
  for l in l1: print l
  for lnr, l in enumerate(l2): 
    if ligand is None or lnr+1 == ligand: 
      has_ens = True
      ll = l.split()
      if ll[0].find(".") > -1 or len(ll) != 7:
        has_ens = False
      if ligand is not None and not has_ens:
        raise ValueError("Ligand %d does not have an ensemble conformation" % ligand)              
      if not has_ens: 
        print l
      else:
        print "    ",
        for f in ll[1:]: print f,
        print
    else:
      print l
