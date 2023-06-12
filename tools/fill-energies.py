import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
ene = open(sys.argv[2])

for h in header: print(h)
stnr = 0
get_next = True
for s in structures:
  while 1:
    if get_next:
      try:
        ll = next(ene)
      except StopIteration:
        raise AssertionError("More structures than energies: %d" % stnr)
    p = ll.find("Energy:")
    if p > -1:
      energy = float(ll[p+len("Energy:"):])
      energy2 = ""
      try:
	      ll = next(ene)
      except StopIteration:
	      pass
      else:
        if ll.find("Energy:") > -1:
          get_next = False
        else:  
          energy2 = ll.rstrip("\n")
      break  
  
  stnr += 1
  l1,l2 = s
  print("#"+str(stnr))
  found = False
  
  for lnr,l in enumerate(l1): 
    if l.startswith("## Energy:"):
      l1[lnr] = "## Energy: " + str(energy)
      found = True
      if len(energy2) and lnr < len(l1) - 1:
        l1[lnr+1] = "## " + energy2
      break    
  if not found: 
    l1.append("## Energy: " + str(energy))
  for l in l1: print(l)
  for l in l2: print(l)

while 1:  
  try:
    next(ene)
  except StopIteration:
    break
  p = ll.find("Energy:")
  if p == -1: continue
  raise AssertionError("More energies than structures: %d" % stnr)
