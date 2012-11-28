#builds an all-atom mode file from a reduced-atom (or CA) mode file
#arguments: mode-aa.py <mode file> <reduced PDB> <all atom PDB>

#the actual numbering may differ between the files
# as long as the residues are in the same order

from __future__ import print_function

import sys
modef, pdbr, pdbaa = sys.argv[1:]

def pdb_to_res(pdbf):
  ret = []
  curr_resid = None
  count = 0
  for l in open(pdbf).readlines():
    if not l.startswith("ATOM"): continue
    resid = l[21:26]
    if resid != curr_resid:
      curr_resid = resid 
      if count: ret.append(count)
      count = 1 
    else:
      count += 1  
  if count > 0: ret.append(count)
  return ret
  
r1 = pdb_to_res(pdbr)
atomsr = sum(r1)
r2 = pdb_to_res(pdbaa)
if len(r1) != len(r2):
  raise Exception("PDBs don't have the same number of residues: %d vs %d" % (len(r1), len(r2)))

def read_modes(modef):
  ret = []
  for l in open(modef).readlines():
    l = l.rstrip("\n")
    ll = l.split()
    if len(ll) == 2 and ll[0] == str(len(ret)+1):
      ret.append([l, []])
      continue
    assert len(ret) > 0  
    ret[-1][1] += ll

  ret0 = ret
  ret = []
  for head, mode in ret0:
    m = []
    for n in range(0,len(mode),3):
      mm = mode[n:n+3]
      m.append(mm)
    ret.append((head,m))
  return ret    
  
modes = read_modes(modef)
for head, mode in modes:
  if len(mode) != atomsr:
    raise Exception("Number of atoms in the mode is not the same as in the PDB: %d vs %d" % (len(mode), atomsr))
  print(head)
  pos = 0
  for count1, count2 in zip(r1, r2):
    vec = mode[pos]
    pr = []
    for n in range(count2):
      pr += vec
    print(" ".join(pr))  
    pos += count1
  

