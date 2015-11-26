"""
Calculate interface RMSD from a multi-model PDB
usage: python irmsd.py 
 [--allatoms] [--allresidues] [--output <file>] [--cutoff <interface cutoff, default 10>]

--allatoms: use all atoms rather than backbone atoms
--allresidues: use also the residues outside the X A interface region 

"""

thresh = 10.0

import sys
import numpy
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib

opt_allatoms = False
opt_allresidues = False

anr = 0
output = None
atomnames = ("CA","C","O","N")

while 1:
  anr += 1
      
  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]

  if arg == "--allatoms": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    opt_allatoms = True
    anr -= 1
    continue
  
  if arg == "--allresidues": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    opt_allresidues = True
    anr -= 1
    continue  
    
  if arg == "--ca": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    atomnames = ("CA",)
    anr -= 1
    continue

  if arg == "--p": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    atomnames = ("P",)
    anr -= 1
    continue
    
  if arg == "--cap": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    atomnames = ("P","CA",)
    anr -= 1
    continue  
       
  if anr <= len(sys.argv)-2 and arg == "--output":
    output = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

  if anr <= len(sys.argv)-2 and arg == "--cutoff":
    thresh = float(sys.argv[anr+1])
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

  if arg.startswith("--"): raise Exception("Unknown option '%s'" % arg)
    
threshsq = thresh * thresh

if len(sys.argv) < 3:
  raise Exception("Please supply a multi-model PDB file and a number of bound PDB files")

boundfiles = sys.argv[2:]

if len(boundfiles) == 1 and opt_allresidues == False:
  raise Exception("Cannot determine the interface for a single PDB")

bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
for p in bounds: 
  p.remove_hydrogens()

allboundatoms = []
for p in bounds:
  for c in p.coordinates():
    allboundatoms.append(c)
allboundatoms = numpy.array(allboundatoms)

selmask = numpy.array([True] * len(allboundatoms))
if not opt_allresidues:
  imask = rmsdlib.build_interfacemask(bounds, thresh)
  selmask = (selmask & imask)
if not opt_allatoms:
  amask = rmsdlib.build_atommask(bounds, atomnames)
  selmask = (selmask & amask)

fboundatoms = allboundatoms[selmask]

nstruc = 0
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')

for unbounds in rmsdlib.read_multi_pdb(sys.argv[1]):
  coor = []  
  for p in unbounds: 
    p.remove_hydrogens()
    for c in p.coordinates():
      coor.append(c)
  coor = numpy.array(coor)
  
  if len(coor) == 0:
    nstruc += 1
    f1.write(str(nstruc)+" None\n")
    continue
  
  rmsdlib.check_pdbs(unbounds, bounds)
  nstruc += 1
  fcoor = numpy.compress(selmask, coor, axis=0)
  irmsd = rmsdlib.fit(fboundatoms,fcoor)[2]
  f1.write(str(nstruc)+" %.3f\n" % irmsd)
