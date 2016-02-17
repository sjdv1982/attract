"""
Calculate interface RMSD per protein partner
Fits the complex on the interface of the reference and then calculates deviations at the interface for each protein partner separately
usage: python irmsd-perligand.py <DAT file> \
 <unbound PDB 1> <bound PDB 1> [<unbound PDB 2> <bound PDB 2>] [...]
 [--allatoms] [--allresidues] [--cutoff <dist cutoff for interface, in A> ]

--allatoms: use all atoms rather than backbone atoms
--allresidues: use also the residues outside the 10 A interface region 

"""
thresh = 10.0

import sys
import numpy
import collectlibpy as collectlib
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib
import math

ensfiles = []
modefile = None
imodefile = None
name = None
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
    
  if anr <= len(sys.argv)-3 and arg == "--ens":
    ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
    sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
    anr -= 3
    continue

  if anr <= len(sys.argv)-2 and arg == "--cutoff":
    thresh = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

  if anr <= len(sys.argv)-2 and arg == "--modes":
    modefile = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
  
  if anr <= len(sys.argv)-2 and arg == "--imodes":
    imodefile = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

# Support for direct IRMSD with iATTRACT files
  if anr <= len(sys.argv)-2 and arg == "--name":
    name = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
  
  if anr <= len(sys.argv)-2 and arg == "--thresh":
    thresh = float(sys.argv[anr+1])
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue 

  if anr <= len(sys.argv)-2 and arg == "--output":
    output = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
  if arg.startswith("--"): raise Exception("Unknown option '%s'" % arg)
    

if len(sys.argv) < 4 or len(sys.argv) % 2:
  raise Exception("Please supply an even number of PDB files (unbound, bound)")

unboundfiles = []
boundfiles = []
for n in range(2, len(sys.argv), 2):
  unboundfiles.append(sys.argv[n])
  boundfiles.append(sys.argv[n+1])

if len(boundfiles) == 1 and opt_allresidues == False:
  raise Exception("Cannot determine the interface for a single PDB")

bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
unbounds = [rmsdlib.read_pdb(f) for f in unboundfiles]

initargs = [sys.argv[1]] + unboundfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
unboundsizes = [len(list(p.atoms())) for p in unbounds]
collectlib.check_sizes(unboundsizes, unboundfiles)

unbound_hmask = rmsdlib.build_hydrogenmask(unbounds)
for p in unbounds + bounds: 
  p.remove_hydrogens()

rmsdlib.check_pdbs(unbounds, bounds)

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

ligandselmask = []
currentligand = 0
currentsel = 0
unboundsizes = [len(list(p.atoms())) for p in unbounds]
for s in unboundsizes:
  count = 0
  for n in selmask[currentligand:currentligand+s]:
    if n:
      count += 1
      
  ligandselmask.append((currentsel,currentsel+count))
  currentsel += count
  currentligand += s
  
fboundatoms = allboundatoms[selmask]

nstruc = 0
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')
while 1:
  sys.stdout.flush()
  if name is not None: 
    newargs = initargs + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
    if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
      break
    collectlib.collect_iattract(newargs)
    
  result = collectlib.collect_next()
  if result: break
  nstruc += 1    
  coor = collectlib.collect_all_coor()
  coor = numpy.compress(unbound_hmask, coor, axis=0)
  fcoor = numpy.compress(selmask, coor, axis=0)
  U = rmsdlib.fit(fcoor,fboundatoms)[0]
  rotfcoor = numpy.array([],ndmin=2)
  L = len(fcoor)
  COM1 = numpy.sum(fcoor,axis=0) / float(L)
  COM2 = numpy.sum(fboundatoms,axis=0) / float(L)
  for c in fcoor:
    rotfcoor = numpy.append(rotfcoor,U.dot(c-COM1))
  
  rotfcoor = numpy.reshape(rotfcoor,(-1,3))
  f1.write(str(nstruc))
  for start, end in ligandselmask:
    sel1 = rotfcoor[start:end]
    sel2 = fboundatoms[start:end]
    sel2 = sel2-COM2
    irmsd = rmsdlib.rmsd(sel1,sel2)
    f1.write(" %.3f" % irmsd)
    
  f1.write("\n")
