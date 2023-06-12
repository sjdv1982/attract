import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
header2,structures2 = read_struc(sys.argv[2])

for h in header: print(h)
for stnr, s1 in enumerate(structures):
  l1,l2 = s1
  try:
    s2 = next(structures2)
  except StopIteration:
    raise AssertionError("Structure file %s contains more structures than %s" % (sys.argv[1], sys.argv[2]))
  l1_orig = s2[0]
     
  print("#"+str(stnr+1))  
  
  for l in l1: 
    print(l)
  for l in l1_orig: 
    if not l.startswith("###"):
      print(l)
  for l in l2: print(l)

try:
  s2 = next(structures2)
except StopIteration:
  pass
else:  
  raise AssertionError("Structure file %s contains more structures than %s" % (sys.argv[2], sys.argv[1]))

