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
import collectlibpy as collectlib
from _read_struc import read_struc
  
def get_selection(atoms, atomnames, use_allatoms):
  
  selatoms = [(n+1,a) for n,a in enumerate(atoms)]
      
  if use_allatoms:
    selatoms = [(n,a) for n,a in selatoms if a[12:16].strip()[0] != "H"]
  else:  
    selatoms = [(n,a) for n,a in selatoms if a[12:16].strip() in atomnames]
  selected = set([n for n,a in selatoms])

  mask = []
  for n in range(len(atoms)):
    mask.append((n+1) in selected)
  return mask

def read_pdb(f):
  ret1 = []
  ret2 = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret1.append((x,y,z))
    ret2.append(l)
  return ret1, ret2
  
  
def rmsd(atoms1, atoms2):
  d = atoms1 - atoms2
  d2 = d * d
  sd = d2.sum(axis=1)
  return numpy.sqrt(sd.mean())

if __name__ == "__main__":  
  import os
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

  unbounds = []
  bounds = []

  
  for n in range(2, len(sys.argv), 2):
    unbounds.append(sys.argv[n])
    bounds.append(sys.argv[n+1])

  initargs = [sys.argv[1], receptor] + unbounds
  if modefile: initargs += ["--modes", modefile]
  if imodefile: initargs += ["--imodes", imodefile]
  ens_receptor = 0
  for nr, ensfile in ensfiles:
    if nr == "1": ens_receptor = 1
    initargs += ["--ens", nr, ensfile]

  collectlib.collect_init(initargs)

  boundatoms = []
  for b in bounds:
    boundatoms.append(read_pdb(b))
  
  unboundatoms = []
  for ub in unbounds:
    unboundatoms.append(read_pdb(ub))
  unboundsizes0 = [len(ub[1]) for ub in unboundatoms]
  
  unboundsizes = []
  start = 0 
  for inr, i in enumerate(collectlib.ieins[:len(unbounds)+1]):
    if inr > 0: unboundsizes.append(i-start)
    else: receptor_offset = i
    start = i
    
  assert unboundsizes0 == unboundsizes, (unboundsizes0, unboundsizes) #ATTRACT and Python disagree about the PDB sizes...
  
  boundmasks = [get_selection(b[1], atomnames, opt_allatoms) for b in boundatoms]
  unboundmasks = [get_selection(ub[1], atomnames, opt_allatoms) for ub in unboundatoms]
  
  
  for bname, ubname, bmask, ubmask in zip(bounds,unbounds,boundmasks,unboundmasks):
    bsize = len([v for v in bmask if v])
    ubsize = len([v for v in ubmask if v])
    if bsize != ubsize:
      raise Exception("Different selected atom numbers: %s: %d, %s: %d" % (ubname, ubsize, bname, bsize))
    if bsize == 0:
      raise Exception("No atoms are selected! %s" % ubname)
 
  allboundatoms = []
  allboundmask = []
  for b, m in zip(boundatoms, boundmasks): 
    allboundatoms.append(numpy.array(b[0]))
    allboundmask.append(numpy.array(m))
  allboundatoms = numpy.concatenate(allboundatoms)  
  allboundmask = numpy.concatenate(allboundmask)
  
  allselbound = allboundatoms[allboundmask]  
  selbound = [numpy.array(at[0])[numpy.array(m)] for at,m in zip(boundatoms, boundmasks)]

  allunboundatoms = []
  allunboundmask = []
  unboundmasks2 = []
  pos = 0
  for ub, m in zip(unboundatoms, unboundmasks): 
    allunboundatoms.append(numpy.array(ub[0]))
    allunboundmask.append(numpy.array(m))
  allunboundatoms = numpy.concatenate(allunboundatoms)  
  allunboundmask = numpy.concatenate(allunboundmask)
  
  for ub, m in zip(unboundatoms, unboundmasks): 
    m2 = [False] * pos + m + [False] * (len(allunboundatoms)-len(m)-pos)
    assert len(m2) == len(allunboundatoms), (len(m2), len(allunboundatoms))
    pos += len(m)
    unboundmasks2.append(numpy.array(m2))
  
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
    
    coor = collectlib.collect_all_coor()
    coor = numpy.array(coor)[receptor_offset:]
    #assert len(coor) == len(allunboundmask)
    allselunbound = coor[allunboundmask]
    f1.write("l-RMSD")
    f1.write(" %.3f" % rmsd(allselbound, allselunbound))
    if len(unboundatoms) > 1:
      for lignr in range(len(unboundatoms)):
        cselunbound = coor[unboundmasks2[lignr]]
        cselbound = selbound[lignr]
        f1.write(" %.3f" % rmsd(cselbound, cselunbound))
    f1.write("\n")
