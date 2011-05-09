import sys
from _read_struc import read_struc
import random
from math import *
header,structures = read_struc(sys.argv[1])

for h in header: print h
stnr = 0
for s in structures:
  stnr += 1
  l1,l2 = s
  for lnr in range(1,len(l2)):
    l = l2[lnr]
    values = [float(v) for v in l.split()]  
    values[0] += 0.35 * (random.random()-0.5)
    values[1] += 0.35 * (random.random()-0.5)/(sin(values[1]+0.1))
    values[2] += 0.35 * (random.random()-0.5)
    for n in 3,4,5:
      values[n] += 7.55 * (random.random()-0.5)
    for n in range(6,len(values)):
      values[n] += 2.5 * (random.random()-0.5) 
    l2[lnr] = "  " + " ".join([("%.6f" % v) for v in values]) 
  print "#"+str(stnr)
  for l in l1: print l
  for l in l2: print l

