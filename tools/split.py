import sys
from _read_struc import read_struc

if len(sys.argv) < 3:
  print >> sys.stderr, "Usage: split.py <DOF file> <file pattern> <nr of splits>"
  sys.exit()

nrsplit = int(sys.argv[3])

pattern = sys.argv[2]
header,structures = read_struc(sys.argv[1])
splitsize = int(float(len(structures))/nrsplit+0.99999)


stnr = 0
totnr = 0
currsplit = 1

f = open("%s-%d" % (pattern,currsplit), "w")
for h in header: print >> f, h
for s in structures:
  totnr += 1
  stnr += 1
  if stnr > splitsize:
    f.close()
    currsplit += 1
    f = open("%s-%d" % (pattern,currsplit), "w")
    for h in header: print >> f, h
    stnr = 1    
  l1,l2 = s
  print >> f, "#"+str(stnr)
  print >> f, "### SPLIT "+str(totnr)
  for l in l1: print >>f, l  
  for l in l2: print >>f, l

