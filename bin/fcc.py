"""
Calculate fraction of common contacts (FCC), see Rodrigues et al. 2012
Does not consider any molecules or parts thereof to be equivalent
usage: python fcc.py <DAT file> \
 <unbound PDB 1> <unbound PDB 2> [<unbound PDB 3>]
 [--cutoff <dist cutoff for interface, in A> ]
"""
thresh = 5.0

import sys
import numpy as np
from scipy.spatial.distance import cdist
import collectlibpy as collectlib
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib

ensfiles = []
modefile = None
imodefile = None
name = None
opt_allatoms = False
opt_allresidues = False

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

# Support for direct FCC with iATTRACT files
  if anr <= len(sys.argv)-2 and arg == "--name":
    name = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue

  if anr <= len(sys.argv)-2 and arg == "--cutoff":
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


unboundfiles = []
for n in range(2, len(sys.argv)):
  unboundfiles.append(sys.argv[n])

if len(unboundfiles) < 2 :
  raise Exception("Cannot determine the contacts for less than two PDBs")

unbounds = [rmsdlib.read_pdb(f) for f in unboundfiles]
residues = [[a.resid for a in u.atoms()] for u in unbounds]

initargs = [sys.argv[1]] + unboundfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
unboundsizes = [len(list(p.atoms())) for p in unbounds]
cum_unboundsizes = [0] + list(np.cumsum(unboundsizes))
collectlib.check_sizes(unboundsizes, unboundfiles)

nstruc = 0
contact_sets = []
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
  coor = collectlib.collect_all_coor()
  nstruc += 1
  current_contact_sets = []
  for n1 in range(len(unbounds)):
    start1, end1 = cum_unboundsizes[n1:n1+2]
    s1 = coor[start1:end1]
    for n2 in range(n1+1, len(unbounds)):
      start2, end2 = cum_unboundsizes[n2:n2+2]
      s2 = coor[start2:end2]
      d = cdist(s1, s2)
      contact1, contact2 = np.where(d<thresh)
      rescontact1 = [residues[n1][a] for a in contact1]
      rescontact2 = [residues[n2][a] for a in contact2]
      rescontacts = set(zip(rescontact1, rescontact2))
      current_contact_sets.append(rescontacts)
  contact_sets.append(current_contact_sets)
for n1 in range(nstruc):
  c1 = contact_sets[n1]
  ncontacts = sum([len(cc) for cc in c1])
  for n2 in range(n1+1, nstruc):
    c2 = contact_sets[n2]
    common = 0
    for i in range(len(c1)):
      cc1, cc2 = c1[i], c2[i]
      common += len(cc1.intersection(cc2))
    print >> f1, n1 + 1, n2 + 1, float(common)/ncontacts
