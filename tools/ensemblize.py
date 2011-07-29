import sys
from math import *
from _read_struc import read_struc
import random

if len(sys.argv) != 5:
  print >> sys.stderr, "Usage: ensemblize <DOF file> <ensemble size> <molecule ID> <mode='all'/'random'>"
  sys.exit()

header,structures = read_struc(sys.argv[1])
ens = int(sys.argv[2])
assert ens > 0

molecule = int(sys.argv[3])
assert molecule > 0

mode = sys.argv[4]
if mode not in ("all", "random"):
  print >> sys.stderr, "Usage: ensemblize <DOF file> <ensemble size> <molecule ID> <mode='all'/'random'>"
  sys.exit()

stnr = 0
for h in header: print h
for s in structures: 
  if mode == "all": ensemble = range(1,ens+1)
  elif mode == "random": ensemble = [random.choice(range(1,ens+1))]
  for e in ensemble:
    stnr += 1
    print "#"+str(stnr)
    l1,l2 = s
    assert len(l2) >= molecule, "Number of molecules is less than %d" % molecule
    head = []
    for l in l2: head.append("")
    head[molecule-1] = str(e) + " "
    for l in l1: print l
    for l,h in zip(l2,head): print h + l
  
