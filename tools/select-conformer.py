#syntax: select-conformer.py <ligand> <conformer>
import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
ligand = int(sys.argv[2])
assert ligand > 0
conformer = int(sys.argv[3])
assert conformer > 0

stnr = 0
for h in header: print h
for s in structures:
  l1,l2 = s
  assert len(l2) >= ligand
  lig = l2[ligand-1]
  ligfields = lig.split()
  if len(ligfields) != 7:
    raise ValueError("Ligand %d, expected 7 fields, read %d" % (ligand, len(ligfields)))
  conf = int(ligfields[0])
  if conf != conformer: continue
  stnr += 1
  try:
    print "#"+str(stnr)
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

