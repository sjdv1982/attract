import sys, glob
from _read_struc import read_struc

if len(sys.argv) < 2 or len(sys.argv) > 3:
  print >> sys.stderr, "Usage: join2.py <file pattern> [--score]"
  sys.exit()
score = False
if len(sys.argv) == 3:
  if sys.argv[2] != "--score":
    print >> sys.stderr, "Usage: join2.py <file pattern> [--score]"
    sys.exit()
  score = True

files = glob.glob(sys.argv[1]+"-*")
nrsplit = 0
while 1:
  fnam = "%s-%d" % (sys.argv[1], nrsplit+1)
  if fnam not in files: break
  nrsplit += 1

if nrsplit == 0: 
  print >> sys.stderr, "Pattern not found"
  sys.exit()

if score:
  for n in range(nrsplit):
    fnam = "%s-%d" % (sys.argv[1], n+1)
    for l in open(fnam).readlines(): 
      print l,
  sys.exit()

allstructures = {}
maxstruc = 0
stnr = 0  
for n in range(nrsplit):
  fnam = "%s-%d" % (sys.argv[1], n+1)
  header0,structures = read_struc(fnam)
  if n == 0: header = header0
  currstruc = None
  currstruc_false = False
  for s in structures:
    stnr += 1
    l1,l2 = s
    print "#"+str(stnr)
    skipline = -1
    for lnr, l in enumerate(l1): 
      if l.startswith("### SPLIT "): 
        skipline = lnr
        break
    for lnr, l in enumerate(l1): 
      if lnr == skipline: continue
      print l
    for l in l2: print l
    
