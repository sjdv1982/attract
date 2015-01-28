import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
selstruc = {}
if sys.argv[2] == "-f":
  selected = [int(l.split()[0]) for l in open(sys.argv[3]).readlines()]
else:  
  selected = [int(v) for v in sys.argv[2:]]
selected = set(selected)

for snr,s in enumerate(structures):
  if snr+1 not in selected: continue
  selstruc[snr+1] =  s

for h in header: print h
for stnr in range(len(selected)):
  snr = selected[stnr]
  l1,l2 = selstruc[snr]
  print "#"+str(stnr+1)
  print "##%d => select" % snr
  try:
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

