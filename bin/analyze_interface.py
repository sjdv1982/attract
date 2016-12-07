"""
Calculate the most 
usage: python stepscore.py <DAT file> <cutoff> <unbound PDBs>
"""

import sys
import numpy as np
import scipy.spatial
from scipy.spatial import cKDTree as KDTree
try:
  KDTree.query_ball_tree
except AttributeError:
  from scipy.spatial import KDTree
  
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

cutoff = float(sys.argv[2])
pdbfiles = sys.argv[3:]
pdbs = [rmsdlib.read_pdb(f) for f in pdbfiles]

initargs = [sys.argv[1]] + pdbfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
pdbsizes = [len(list(p.atoms())) for p in pdbs]
collectlib.check_sizes(pdbsizes, pdbfiles)
resids = []
for pdb in pdbs:
  cresids = []
  for residue in pdb.residues():
    resid = residue[0].resid
    for n in range(len(residue)):
      cresids.append(resid)
  resids.append(cresids)


f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')
  
nstruc = 0
rescounts = [{} for n in range(len(pdbs))]  
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
  for n1 in range(len(pdbs)):
    c1 = coors[n1]
    tree1 = KDTree(c1)
    resids1 = resids[n1]
    rescounts1 = rescounts[n1]
    for n2 in range(n1+1, len(pdbs)):
      c2 = coors[n2]
      tree2 = KDTree(c2)
      resids2 = resids[n2]
      rescounts2 = rescounts[n2]
      pairs = tree1.query_ball_tree(tree2, cutoff)
      for p1,pp in enumerate(pairs):
        if not len(pp): continue
        resid1 = resids1[p1]
        if resid1 not in rescounts1:
          rescounts1[resid1] = 0
        rescounts1[resid1] += len(pp)
        for p2 in pp:
          resid2 = resids2[p2]
          if resid2 not in rescounts2:
            rescounts2[resid2] = 0
          rescounts2[resid2] += 1
avg_contacts = [sum(c.values())/float(nstruc) for c in rescounts]
#from pprint import pprint
#pprint(rescounts) 
#print(sum(rescounts[0].values()), sum(rescounts[1].values()), avg_contacts)

print >> f1, "Partner\tresid\tcontacts"
for n in range(len(pdbs)):
  rescount = rescounts[n]
  for r in rescount:
      rescount[r] = float(rescount[r])/nstruc
  res = list(rescount.keys())
  res.sort(key=lambda k:rescount[k],reverse=True)
  accum = 0  
  last = None
  for rnr, r in enumerate(res):
      c = rescount[r]
      if c != last:
        tie = 1
        for r2 in res[rnr+1:]: 
           if rescount[r2] != c:
             break
           tie += 1
        accum += c * tie
        if accum >= avg_contacts[n]:
          break      
        last = c
      if r[0] == " ":
          r = "_" + r[1:]
      print >> f1, n+1, r,  c      
