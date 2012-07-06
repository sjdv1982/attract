import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
nrclones = int(sys.argv[2])

energies = []
stnr = 0
stnr2 = 0
for h in header: print h
for l1,l2 in structures:
  stnr += 1
  for n in range(nrclones):
    stnr2 += 1
    print "#"+str(stnr2)
    print "##"+str(stnr) + " => clone"
    for l in l1: print l
    for l in l2: print l
