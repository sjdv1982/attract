import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])


energies = []
for ll in open(sys.argv[2]).readlines():
  p = ll.find("Energy:")
  if p > -1:
    e = float(ll[p+len("Energy:"):])
    energies.append(e)
assert len(energies) == len(structures), (len(energies), len(structures))
strucs = zip(range(1,len(structures)+1), energies, structures)

for h in header: print h
stnr = 0
for st in strucs:
  r,e,s = st
  stnr += 1
  l1,l2 = s
  print "#"+str(stnr)
  found = False
  for l in l1: 
    if l.startswith("## Energy:"):
      print "## Energy:", e
      found = True
      continue
    print l
  if not found: print "## Energy:", e
  for l in l2: print l

