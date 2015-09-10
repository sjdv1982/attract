"""
Sets all pivots of a .dat file to zero
"""
import sys
from _read_struc import read_struc
header,strucs = read_struc(sys.argv[1])
from euler2rotmat import euler2rotmat

from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option("--ens", 
                  action="append", dest="ens", type="int")
options, args = parser.parse_args()
ens = options.ens
if ens is None: ens = []

pivots = []
for h in header: 
  if not h.startswith("#pivot"):
    h = h.rstrip()
    if h.startswith("#centered"): assert h.endswith(" false"), h
    print h
    continue
  assert not h.startswith("#pivot auto"), h
  hh = h.split()
  assert hh[1] == str(len(pivots)+1), h
  assert len(hh) == 5, h
  pivot = [float(v) for v in hh[2:5]]
  pivots.append(np.array(pivot))
  print "#pivot %d 0 0 0" % len(pivots)
  
stnr = 0
for st in strucs:
  stnr += 1
  l1,l2 = st
  print "#"+str(stnr)
  try:
    for l in l1: print l
    for lnr, l in enumerate(l2): 
      ll = l.split()
      has_ens = ((lnr+1) in ens)
      assert len(ll) == 6 + has_ens
      if has_ens:
        ensc = ll[0]
        ll = ll[1:]
      values = [float(vv) for vv in ll]
      rotmat = np.array(euler2rotmat(*values[:3]))
      p = pivots[lnr]
      pp = (-p * rotmat).sum(axis=1) + p
      trans = np.array(values[3:])
      trans += pp
      if has_ens:
        print ensc,        
      print values[0], values[1], values[2], trans[0], trans[1], trans[2]
  except IOError:  
    sys.exit()

