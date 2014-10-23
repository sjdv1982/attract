"""
Adds TER statements to a PDB by detecting chain breaks:
- For protein: C-N distance of >1.5 A
- For DNA/RNA: P-O3'distance of >2.0 A
TER statements are also added between protein and DNA/RNA
Existing TER statements and chain IDs are ignored

Author: Sjoerd de Vries, Technische Universitaet Muenchen
"""

import sys, os
from math import *


filename = sys.argv[1]
data = open(filename).readlines()
data = [l for l in data if l.startswith("ATOM") ]
residues = []

cres = None
for l in data:
  resid = l[21:26]
  if resid != cres:
    residues.append([])
    cres = resid
  residues[-1].append(l)

def parse_res(r):
  ret = {}
  for l in r:
    name = l[12:16].strip()
    x = float(l[30:38])
    y = float(l[38:46])
    z = float(l[46:54])
    ret[name] = x,y,z
  return ret

lastpos = None
lastname = None
for res in residues:
  atoms = parse_res(res)
  name = None
  nextname = None
  if "N" in atoms and "C" in atoms:
    name = "N"
    nextname = "C"
  elif "O3'" in atoms and "P" in atoms:  
    name = "P"
    nextname = "O3'"
  has_ter = False
  if lastname is not None:
    has_ter = True
    if lastname is not None and name is not None:
      has_ter = False
      pos = atoms[name]
      dx = lastpos[0]-pos[0]
      dy = lastpos[1]-pos[1]
      dz = lastpos[2]-pos[2]
      dissq = dx**2+dy**2+dz**2
      if (lastname, name) == ("C", "N"):
        lim = 1.5
      elif (lastname, name) == ("O3'", "P"):
        lim = 2.0
      else:
        lim = 0
      if dissq > lim*lim: has_ter = True
  if has_ter: print "TER"
  for l in res: print l,
  
  lastname = nextname
  lastpos = None
  if lastname is not None:
    lastpos = atoms[lastname]
  print "TER"  