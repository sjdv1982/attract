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

options, args = parser.parse_args()

start = 1
if options.fix_receptor: start = 0
dori = options.ori
dtrans = options.trans
dmode = options.mode

random.seed(options.seed)

import sys
from _read_struc import read_struc
import random
from math import *
header,structures = read_struc(args[0])

for h in header: print h
stnr = 0
for s in structures:
  stnr += 1
  l1,l2 = s
  for lnr in range(start,len(l2)):
    l = l2[lnr]
    values = [float(v) for v in l.split()]  
    values[0] += dori * (random.random()-0.5)
    values[1] += dori * (random.random()-0.5)/(sin(values[1]+0.1))
    values[2] += dori * (random.random()-0.5)
    for n in 3,4,5:
      values[n] += dtrans * (random.random()-0.5)
    for n in range(6,len(values)):
      values[n] += dmode * (random.random()-0.5) 
    l2[lnr] = "  " + " ".join([("%.6f" % v) for v in values]) 
  print "#"+str(stnr)
  for l in l1: print l
  for l in l2: print l

