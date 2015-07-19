import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
selected = [int(v) for v in sys.argv[2:]]

pivots = []
header2 = []
for h in header:
  hh = h.split()
  if h.startswith("##"):
    print h
    continue
  if h.startswith("#pivot"):
    assert not hh[1] == "auto"
    assert hh[1] == str(len(pivots)+1)
    assert len(hh) == 5
    pivots.append(hh[2:5])    
  else:
    header2.append(h)

count = 0
for pnr in range(len(pivots)):
  if pnr+1 not in selected: continue
  count += 1
  print "#pivot", count, 
  for pv in pivots[pnr]: print pv,
  print
for h in header2: print h

for snr,s in enumerate(structures):
  l1,l2 = s
  print "#"+str(snr+1)
  try:
    for l in l1: print l
    for lnr, l in enumerate(l2): 
      if lnr+1 in selected:
        print l
  except IOError:  
    sys.exit()

