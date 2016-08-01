import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
rev = False
if len(sys.argv) > 2:
  assert len(sys.argv) == 3
  assert sys.argv[2].startswith("--rev")
  rev = True

energies = []
for l1,l2 in structures:
  for ll in l1:
    if ll.startswith("## Energy:"):
      ee = ll[10:].strip()
      if ee.startswith("nan"):
        e = 99999999999999
      else:
        e = float(ee)
      energies.append(e)
assert len(energies) == len(structures)
strucs = zip(range(1,len(structures)+1), energies, structures)

strucs.sort(key=lambda v:v[1])
if rev: strucs.reverse()
for h in header: print h
stnr = 0
for st in strucs:
  r,e,s = st
  stnr += 1
  l1,l2 = s
  try:  
    print "#"+str(stnr)
    print "## "+str(r) + " => sort"
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

