import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
energy_threshold = float(sys.argv[2])
structures = list(structures)
rev = False
if len(sys.argv) > 3:
  assert len(sys.argv) == 4
  assert sys.argv[3].startswith("--rev")
  rev = True

energies = []
for l1,l2 in structures:
  for ll in l1:
    if ll.startswith("## Energy:"):
      e = float(ll[10:])
      energies.append(e)
assert len(energies) == len(structures)
strucs = zip(range(1,len(structures)+1), energies, structures)

for h in header: print h
stnr = 0
for st in strucs:
  r,e,s = st  
  l1,l2 = s
  if rev:
    if e < energy_threshold: continue
  else:
    if e > energy_threshold: continue
  stnr += 1    
  print "#"+str(stnr)
  print "##"+str(r) + " => filter-energy"
  try:
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

