import sys

if len(sys.argv) not in (9,10):
  print >> sys.stderr, "Usage: air.py <protein 1 active residue list> <protein 1 passive residue list>"  
  print >> sys.stderr, "              <protein 1 reduced PDB> <protein 1 residue mapping file>"  
  print >> sys.stderr, "              <protein 2 active residue list> <protein 2 passive residue list>"  
  print >> sys.stderr, "              <protein 2 reduced PDB> <protein 2 residue mapping file>"    
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
    ret[ll[0]] = int(ll[1])
  return ret
    
act1 = [l.strip() for l in open(sys.argv[1]).readlines()]
pass1 = [l.strip() for l in open(sys.argv[2]).readlines()]
pdb1,names1,pdb1len = load_pdb(sys.argv[3])
map1 = load_map(sys.argv[4])
act1a = [map1[r] for r in act1]
pass1a = [map1[r] for r in pass1]

act2 = [l.strip() for l in open(sys.argv[5]).readlines()]
pass2 = [l.strip() for l in open(sys.argv[6]).readlines()]
pdb2,names2,pdb2len = load_pdb(sys.argv[7])
map2 = load_map(sys.argv[8])
act2a = [map2[r] for r in act2]
pass2a = [map2[r] for r in pass2]

actpass1 = []
for a in act1a+pass1a:
  actpass1 += pdb1[a]
print "A_actpass", len(actpass1),
for nr in actpass1: print nr,
print

actpass2 = []
for a in act2a+pass2a:
  actpass2 += pdb2[a]
print "B_actpass", len(actpass2),
for nr in actpass2: print nr+pdb1len,
print
  
for r in act1:
  nrs = pdb1[map1[r]]
  print "A_"+names1[map1[r]]+r,len(nrs),
  for nr in nrs: print nr,
  print
for r in act2:
  nrs = pdb2[map2[r]]
  print "B_"+names2[map2[r]]+r,len(nrs),
  for nr in nrs: print nr+pdb1len,
  print

print

noecv = 0.5
if len(sys.argv) == 10: noecv = float(sys.argv[9])
params = 2, 2.0, 1.0, 2.0, noecv


for r in act1:
  print "A_"+names1[map1[r]]+r, "B_actpass", 
  for p in params: print p,
  print
for r in act2:
  print "B_"+names2[map2[r]]+r, "A_actpass", 
  for p in params: print p,
  print
  
