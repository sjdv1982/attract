"""
Generates HADDOCK/CNS .tbl file from a list of HADDOCK-style active and passive residues
Uses the generate_tbl code from the HADDOCK server
Author: Sjoerd de Vries
"""
def generate_tbl(act_1, pass_1, act_2, pass_2, segid1, segid2, dist): 
  ret = ""
  ret += "! HADDOCK AIR restraints for 1st partner\n"
  actpass2 = act_2 + pass_2
  if len(actpass2) > 0:
    for r1 in act_1:
      ret += "!\n"
      ret += "assign ( resid %s  and segid %s)\n" % (r1, segid1)
      ret += "       (\n"
      first = 1
      for r2 in actpass2:
        if (first == 0):
          ret += "     or\n"
        else: first = 0
        ret += "        ( resid %s  and segid %s)\n" % (r2, segid2)
      ret += "       )  2.0 2.0 0.0\n"
  
  ret += "!\n"
  ret += "! HADDOCK AIR restraints for 2nd partner\n"
  actpass1 = act_1 + pass_1
  if len(actpass1) > 0:
    for r1 in act_2:
      ret += "!\n"
      ret += "assign ( resid %s  and segid %s)\n" % (r1, segid2)
      ret += "       (\n"
      first = 1
      for r2 in actpass1:
        if (first == 0):
          ret += "     or\n"
        else: first = 0
        ret += "        ( resid %s  and segid %s)\n" % (r2, segid1)
      ret += "       )  %f %f 0.0\n" % (dist, dist)
  return ret

import sys

if len(sys.argv) not in (7,8):
  print >> sys.stderr, "Usage: generate_tbl.py <protein 1 active residue list> <protein 1 passive residue list>"  
  print >> sys.stderr, "              <protein 2 active residue list> <protein 2 passive residue list>"  
  print >> sys.stderr, "              <chain 1> <chain 2> [minimum AIR distance]" 
  sys.exit()

def load_pdb(pdb):
  ret = [{}, {}, 0]
  curr = None
  for lnr,l in enumerate(open(pdb).readlines()):
    if not l.startswith("ATOM "): continue
    ret[2] += 1
    aa = int(l[22:27])
    if aa is not curr:
      curr = aa
      ret[0][curr] = []
      ret[1][curr] = l[17]+l[18:20].lower()
    ret[0][curr].append(lnr+1)
  return ret
    
def load_map(mapfile):
  ret = {}
  for l in open(mapfile).readlines():
    if not len(l.strip()): continue
    ll = l.split()
    try:
      ret[ll[0]] = int(ll[1])
    except:
      raise ValueError(l)
  return ret
    
act1 = [l.strip() for l in open(sys.argv[1]).readlines()]
pass1 = [l.strip() for l in open(sys.argv[2]).readlines()]

act2 = [l.strip() for l in open(sys.argv[3]).readlines()]
pass2 = [l.strip() for l in open(sys.argv[4]).readlines()]

chain1 = sys.argv[5]
chain2 = sys.argv[6]

dist = 2.0 
if len(sys.argv) >= 8: dist = float(sys.argv[7])
tbl = generate_tbl(act1, pass1, act2, pass2, chain1, chain2, dist)
print tbl

