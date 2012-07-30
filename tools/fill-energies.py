import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])


energies = []
for ll in open(sys.argv[2]).readlines():
  p = ll.find("Energy:")
  if p > -1:
    e = float(ll[p+len("Energy:"):])
    energies.append(e)

for h in header: print h
stnr = 0
for e in energies:
  try:
    s = structures.next()
  except StopIteration:
    raise IndexError(stnr)
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

try:
  structures.next()
except StopIteration:
  pass  
else:
  raise AssertionError("More structures than energies: %d" % stnr)
