import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
selstruc = int(sys.argv[2])
rev = False
if len(sys.argv) > 3:
  assert len(sys.argv) == 4
  assert sys.argv[3].startswith("--rev")
  rev = True

energies = []
combi = []
for l1,l2 in structures:
  for l in l1:    
    if l.startswith("## Energy:"):
      ee = l[10:].strip()
      if ee.startswith("nan"):
        e = 99999999999999
      else:
        e = float(ee)
      energies.append(e)
      continue
    ll = l.split()
    if ll[0] == "##" and ll[1].find("combine") > -1:
      ccombi = [int(v) for v in ll[2:]]
      combi.append(ccombi)
assert len(energies) == len(structures)

origin = {}
selected = set()
for stnr, c in enumerate(combi):
  cc = set(c)
  if len(cc) == 1: #original structure
    selected.add(stnr)
    continue
  for v in cc:
    if v not in origin: origin[v] = []
    origin[v].append(stnr)

if selstruc < len(structures):  
  allowed = {}
  sorigin = sum([len(v) for v in origin.values()])
  csorigin = 0
  callowed = 0
  oristrucs = list(reversed(sorted(origin.keys())))
  selstruc_remaining = selstruc - len(selected)
  for o in oristrucs:
    csorigin += len(origin[o])
    callowed_new = int(float(csorigin)/sorigin * selstruc_remaining+0.5)
    allowed[o] = callowed_new - callowed
    callowed = callowed_new
    
  for o in oristrucs:
    candidates = [(stnr, energies[stnr]) for stnr in origin[o] if stnr not in selected]
    candidates.sort(key=lambda c: c[1])
    if rev: candidates.reverse()
    for c in candidates[:allowed[o]]:
      selected.add(c[0])
else:
  selected = set(range(len(structures)))

strucs = zip(range(1,len(structures)+1), energies, structures)
strucs.sort(key=lambda v:v[1])
if rev: strucs.reverse()
for h in header: print h
stnr = 0
for st in strucs:
  r,e,s = st
  if r-1 not in selected: continue
  stnr += 1
  l1,l2 = s
  print "#"+str(stnr)
  print "##"+str(r) + " => topcombine"
  try:
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()
