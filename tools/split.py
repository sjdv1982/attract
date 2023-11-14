from __future__ import print_function
import sys
from _read_struc import read_struc

if len(sys.argv) < 3:
  print("Usage: split.py <DOF file> <file pattern> <nr of splits>", file=sys.stderr)
  sys.exit()

datfile = sys.argv[1]
nrsplit = int(sys.argv[3])

pattern = sys.argv[2]

#determine number of structures
for l in reversed(open(datfile).readlines()):
  if not l.startswith("#"): continue
  if l.startswith("##"): continue
  nrstruc = int(l[1:])
  break
header,structures = read_struc(datfile)
splitsize = int(float(nrstruc)/nrsplit+0.99999)

stnr = 0
totnr = 0
currsplit = 1

filename = "%s-%d" % (pattern,currsplit)
print(filename)
f = open(filename, "w")
for h in header: print(h, file=f)
for s in structures:
  totnr += 1
  stnr += 1
  if stnr > splitsize:
    f.close()
    currsplit += 1
    filename = "%s-%d" % (pattern,currsplit)
    print(filename)
    f = open(filename, "w")
    for h in header: print(h, file=f)
    stnr = 1
  l1,l2 = s
  print("#"+str(stnr), file=f)
  print("### SPLIT "+str(totnr), file=f)
  for l in l1: print(l, file=f)
  for l in l2: print(l, file=f)
