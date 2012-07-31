import sys
from math import *
from _read_struc import read_struc
import random

if len(sys.argv) not in (7,8):
  print >> sys.stderr, "Usage: locrest <DOF file> <molecule ID> <x y z> <margin> [--has-ensemble]"
  sys.exit()

has_ensemble = False
if len(sys.argv) == 8:
  assert sys.argv[7] == "--has-ensemble"
  has_ensemble = True

pos = 6
if has_ensemble: pos = 7

header,structures = read_struc(sys.argv[1])

molecule = int(sys.argv[2])
assert molecule > 0

loc = [float(a) for a in sys.argv[3:6]]
margin = float(sys.argv[6])

stnr = 0
for h in header: print h
for s in structures: 
  stnr += 1
  print "#"+str(stnr)
  l1,l2 = s
  assert len(l2) >= molecule, "Number of molecules is less than %d" % molecule
  for l in l1: print l
  for lnr,l in enumerate(l2):
    if lnr != molecule-1:
      print "   " + l.lstrip()
      continue
    ll = l.split()
    assert len(ll) >= pos, (len(ll), pos)
    print "  ",
    for n in range(pos): 
      print ll[n],
    for n in range(3): 
      r = 1-2*random.random()
      target = loc[n] + margin * r
      print target,
    for n in range(pos,len(ll)): 
      print ll[n],
    print      
  
