"""
Calculate ligand RMSD
usage: python lrmsd.py <DAT file> \
 <unbound PDB 1> <bound PDB 1> [<unbound PDB 2> <bound PDB 2>] [...]
 [--allatoms] [--ca] [--p]

--allatoms: use all atoms rather than backbone atoms
--ca: use CA atoms rather than backbone atoms
--p: use P atoms rather than backbone atoms (nucleic acids)
--receptor, --imodes, --modes, --name, --ens, --output: ...
"""
import sys

import numpy
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import collectlibpy as collectlib
from _read_struc import read_struc
import rmsdlib
  
ensfiles = []
modefile = None
imodefile = None
name = None
opt_allatoms = False
atomnames = ("CA", "C", "O", "N")
receptor = "/dev/null"

anr = 0
output = None
while 1:
  anr += 1
      
  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]

  if arg == "--allatoms": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    opt_allatoms = True
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

if arg == "--nucleic-acid": 
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    atomnames = ("P","CA",)
    anr -= 1
    continue  
        
  if anr <= len(sys.argv)-3 and arg == "--ens":
    ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
    sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
    anr -= 3
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

  if anr <= len(sys.argv)-2 and arg == "--receptor":
    receptor = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

# Suppport for direct IRMSD with iATTRACT files
  if anr <= len(sys.argv)-2 and arg == "--name":
    name = sys.argv[anr+1]
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
    
bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
unbounds = [rmsdlib.read_pdb(f) for f in unboundfiles]

initargs = [sys.argv[1], receptor] + unboundfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
ens_receptor = 0
for nr, ensfile in ensfiles:
  if nr == "1": ens_receptor = 1
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
unboundsizes = [len(list(p.atoms())) for p in unbounds]
receptor_offset = collectlib.ieins[0]
collectlib.check_sizes([receptor_offset] + unboundsizes, [receptor] + unboundfiles)


if opt_allatoms:
  unbound_amask = rmsdlib.build_hydrogenmask(unbounds)
else:
  unbound_amask = rmsdlib.build_atommask(unbounds, atomnames)

  
unbound_amasks_ligand = []
for unr, u in enumerate(unbounds):
  if opt_allatoms:
    uu = rmsdlib.build_hydrogenmask([u])
  else:
    uu = rmsdlib.build_atommask([u], atomnames)
  offset = 0
  if unr: offset = sum(unboundsizes[:unr])
  offsetmask = numpy.zeros(offset, dtype="bool")
  uu = numpy.concatenate((offsetmask, uu),axis=0)
  unbound_amasks_ligand.append(uu)
  
for p in unbounds + bounds: 
  if opt_allatoms:
    p.remove_hydrogens()
  else:
    p.select_atoms(atomnames)

rmsdlib.check_pdbs(unbounds, bounds)

allboundatoms = []
boundatoms = []
for p in bounds:
  b = []
  for c in p.coordinates():
    allboundatoms.append(c)
    b.append(c)
  boundatoms.append(numpy.array(b))
allboundatoms = numpy.array(allboundatoms)

nstruc = 0
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')
h, strucs = read_struc(sys.argv[1])
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
  
  l1, l2 = strucs.next()
  ll = [float(v) for v in l2[0].split()[ens_receptor:ens_receptor+6]]
  for v in ll:
    if abs(v)> 0.001:
      raise ValueError("Structures have not yet been fitted")
  
  f1.write("l-RMSD")
  coor = collectlib.collect_all_coor()
  coor = numpy.array(coor)[receptor_offset:]
  fcoor = numpy.compress(unbound_amask, coor, axis=0)
  rmsd = rmsdlib.rmsd(allboundatoms,fcoor)
  f1.write(" %.3f" % rmsd)
  if len(unbounds) > 1:
    for n in range(len(unbounds)):
      fcoor = numpy.compress(unbound_amasks_ligand[n], coor, axis=0)      
      rmsd = rmsdlib.rmsd(boundatoms[n],fcoor)
      f1.write(" %.3f" % rmsd)
  f1.write("\n")
