"""
Calculate Alex step score
usage: python stepscore.py <DAT file> <unbound PDBs>
"""

import sys
import numpy as np
import scipy.spatial
from scipy.spatial import cKDTree as KDTree
import collectlibpy as collectlib
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib

ensfiles = []
modefile = None
imodefile = None
name = None

anr = 0
output = None
while 1:
  anr += 1

  if anr > len(sys.argv)-1: break
  arg = sys.argv[anr]

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

# Support for direct IRMSD with iATTRACT files
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

pdbfiles = sys.argv[2:]
pdbs = [rmsdlib.read_pdb(f) for f in pdbfiles]

initargs = [sys.argv[1]] + pdbfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
pdbsizes = [len(list(p.atoms())) for p in pdbs]
collectlib.check_sizes(pdbsizes, pdbfiles)
atomtypes = [[int(a.line[57:59]) for a in p.atoms()] for p in pdbs]

params = np.loadtxt(os.environ["ATTRACTTOOLS"]+"/../MC_gaa_step10.par")
energyblocks = []
for n1 in range(len(pdbs)):
  a1 = atomtypes[n1]
  for n2 in range(n1+1, len(pdbs)):
    a2 = atomtypes[n2]
    energyblock = np.zeros((len(a1), len(a2)),dtype=params.dtype)
    for i in range(len(a1)):
      t1 = a1[i]
      for j in range(len(a2)):
        t2 = a2[j]
        energyblock[i,j] = params[t1-1, t2-1]
    energyblocks.append(energyblock)


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
  pdbsizes2 = np.cumsum([0] + pdbsizes)
  coors = [coor[pdbsizes2[n]:pdbsizes2[n+1]] for n in range(len(pdbs))]
  energy = 0
  eblock = 0
  for n1 in range(len(pdbs)):
    c1 = coors[n1]
    tree1 = KDTree(c1)
    for n2 in range(n1+1, len(pdbs)):
      c2 = coors[n2]
      tree2 = KDTree(c2)
      energyblock = energyblocks[eblock]
      pairs = tree1.query_ball_tree(tree2, 10)
      ene = sum([sum(e[p]) for e,p in zip(energyblock,pairs) if len(p)])
      energy += ene

      eblock += 1


  f1.write("%.3f\n" % energy)
