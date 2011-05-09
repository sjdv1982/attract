import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
header2,structures2 = read_struc(sys.argv[2])

assert len(structures) == len(structures2)

strucs = zip(structures,structures2)

for h in header: print h
stnr = 0
for st in strucs:
  s,s2 = st
  stnr += 1
  l1,l2 = s
  ll1,ll2 = s2
  print "#"+str(stnr)
  for l in l1: 
    if l not in ll1: print l
  for l in ll1: print l
  for l in l2: print l

