import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
header2,structures2 = read_struc(sys.argv[2])

for h in header: print(h)
st2nr = 0
for stnr, s1 in enumerate(structures):
  l1,l2 = s1
  dr = None
  for l in l1:
    ll = l.split()
    if len(ll) == 4 and ll[3] == "deredundant" and ll[0] == "##" and ll[2] == "=>":
      dr = int(ll[1])
      break
  else:
    raise ValueError("Cannot find deredundant marker for structure %d" % (stnr+1))
  while st2nr < dr:
    try:
      s2 = next(structures2)
      st2nr += 1
    except StopIteration:
      raise AssertionError("Cannot find structure: %d => deredundant %d" % (dr, stnr))
  l1_orig = s2[0]
     
  print("#"+str(stnr+1))  
  
  for l in l1: 
    print(l)
  for l in l1_orig: 
    if not l.startswith("###"):
      print(l)
  for l in l2: print(l)
