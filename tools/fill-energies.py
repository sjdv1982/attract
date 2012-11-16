import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])


energies = []
energies2 = []
lines = open(sys.argv[2]).readlines()
for llnr,ll in enumerate(lines):
  p = ll.find("Energy:")
  if p > -1:
    e = float(ll[p+len("Energy:"):])
    energies.append(e)
    try:
      ll2 = lines[llnr+1]
    except IndexError:
      continue
    if len(ll2.split()) == 6: energies2.append(ll2.rstrip("\n"))

for h in header: print h
stnr = 0
for enr, e in enumerate(energies):
  try:
    s = structures.next()
  except StopIteration:
    raise IndexError(stnr)
  stnr += 1
  l1,l2 = s
  print "#"+str(stnr)
  found = False
  
  for lnr,l in enumerate(l1): 
    if l.startswith("## Energy:"):
      l1[lnr] = "## Energy: " + str(e)
      found = True
      if len(energies2) == len(energies) and lnr < len(l1) - 1:
        l1[lnr+1] = "## " + energies2[enr]
      break    
  if not found: l1.append("## Energy: " + str(e))
  for l in l1: print l
  for l in l2: print l

try:
  structures.next()
except StopIteration:
  pass  
else:
  raise AssertionError("More structures than energies: %d" % stnr)
