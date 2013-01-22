import time, random

from optparse import OptionParser

parser = OptionParser()

parser.add_option("--fix-receptor",
                  action="store_true", dest="fix_receptor", default=False,
                  )
parser.add_option("--ori",
                  action="store", dest="ori", default = 3.5, type="float",
                  )
parser.add_option("--trans", "--tra",
                  action="store", dest="trans", default = 7.55, type="float",
                  )
parser.add_option("--mode", 
                  action="store", dest="mode", default = 2.5, type="float",
                  )
parser.add_option("--seed", 
                  action="store", dest="seed", default = time.time(), type="int",
                  )

parser.add_option("--morph", 
                  action="append", dest="morph", type="int",
                  )
parser.add_option("--ens", 
                  action="append", dest="ens", type="int",
                  )
parser.add_option("--clone", 
                  action="store", dest="clone", default = -1, type="int",
                  )
parser.add_option("--keepfirst", 
                  action="store_true", dest="keepfirst", default = False,
                 )

options, args = parser.parse_args()
morph, ens = options.morph, options.ens
if morph is None: morph = []
if ens is None: ens = []
for m in morph: assert m not in ens, m

start = 0
if options.fix_receptor: start = 1
dori = options.ori
dtrans = options.trans
dmode = options.mode
clone = options.clone
keepfirst = options.keepfirst

random.seed(options.seed)

import sys
from _read_struc import read_struc
import random
from math import *
header,structures = read_struc(args[0])

for h in header: print h
stnr = 0
for s in structures:
  clonenr = 0
  while 1:
    stnr += 1
    l1,l2 = s
    l1,l2 = list(l1), list(l2) #copy    
    for lnr in range(start,len(l2)):
      l = l2[lnr]
      values = [float(v) for v in l.split()]  
      is_morph, is_ens = False, False
      if lnr+1 in morph: is_morph = True
      elif lnr+1 in ens: is_ens = True
      if is_morph:
        cmorph = values[0]
        values = values[1:]
      elif is_ens:
        cens = values[0]
        values = values[1:]
      values[0] += dori * (random.random()-0.5)
      #values[1] += dori * (random.random()-0.5)/(sin(values[1]+0.1))
      values[1] += dori * (random.random()-0.5)
      values[2] += dori * (random.random()-0.5)
      for n in 3,4,5:
        values[n] += dtrans * (random.random()-0.5)
      for n in range(6,len(values)):
        values[n] += dmode * (random.random()-0.5) 
      if is_morph:
        cmorph += dmode * (random.random()-0.5) 
        if cmorph < 0: cmorph = 0
        values = [cmorph] + values
      elif is_ens:
        values = [cens] + values
      l2[lnr] = "  " + " ".join([("%.6f" % v) for v in values]) 
    if clonenr == 0 and keepfirst: 
      l1,l2 = s
    print "#"+str(stnr)
    for l in l1: print l
    for l in l2: print l
    if clone == -1: break
    clonenr += 1
    if clonenr >= clone: break
