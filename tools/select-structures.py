import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
selected = [int(v) for v in sys.argv[2:]]


for h in header: print h
stnr = 0
for snr,s in enumerate(structures):
  if snr+1 not in selected: continue
  stnr += 1
  l1,l2 = s
  print "#"+str(stnr)
  print "##%d => select" % (snr+1)
  try:
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

