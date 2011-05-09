import sys
from _read_struc import read_struc
import random
from math import *
header,structures = read_struc(sys.argv[1])
header2,structures2 = read_struc(sys.argv[2])
##temperature = float(argv[3])

def select(e1,e2):
  return e2 < e1

def get_energy(l1):
  for l in l1:
    if l.startswith("## Energy: "):
      energy = float(l[11:])
      return energy

for h in header: print h
stnr = 0
strucs = zip(structures,structures2)

for s1,s2 in strucs:
  stnr += 1
  l1,l2 = s1
  ll1,ll2 = s2
  
  energy1 = get_energy(l1)
  energy2 = get_energy(ll1)

  if select(energy1,energy2):
    l1 = ll1
    l2 = ll2

  print "#"+str(stnr)
  for l in l1: print l
  for l in l2: print l

