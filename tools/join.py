import sys, glob
from _read_struc import read_struc

if len(sys.argv) < 2 or len(sys.argv) > 3:
  print("Usage: join.py <file pattern> [--score]", file=sys.stderr)
  sys.exit()
score = False
if len(sys.argv) == 3:
  if sys.argv[2] != "--score":
    print("Usage: join.py <file pattern> [--score]", file=sys.stderr)
    sys.exit()
  score = True

files = glob.glob(sys.argv[1]+"-*")
nrsplit = 0
while 1:
  fnam = "%s-%d" % (sys.argv[1], nrsplit+1)
  if fnam not in files: break
  nrsplit += 1

if nrsplit == 0: 
  print("Pattern not found", file=sys.stderr)
  sys.exit()

if score:
  for n in range(nrsplit):
    fnam = "%s-%d" % (sys.argv[1], n+1)
    for l in open(fnam).readlines(): 
      assert not l.startswith("#") #must be .score file, not .dat file!!
      print(l, end=' ')
  sys.exit()

allstructures = {}
maxstruc = 0
stnr = 0  
for n in range(nrsplit):
  fnam = "%s-%d" % (sys.argv[1], n+1)
  header0,structures = read_struc(fnam)
  if n == 0: 
    for h in header0:
      print(h)
  currstruc = None
  currstruc_false = False
  for s in structures:
    stnr += 1
    l1,l2 = s
    print("#"+str(stnr))
    skipline = -1
    for lnr, l in enumerate(l1): 
      if l.startswith("### SPLIT "): 
        try:
          currstruc = int(l[len("### SPLIT"):])
        except:
          currstruc = l[len("### SPLIT"):]
          currstruc_false = True
        skipline = lnr
        break

    if currstruc is None:
      print("Structure has no SPLIT number", file=sys.stderr)
      print(fnam, file=sys.stderr)
      print("#"+str(stnr), file=sys.stderr)
      for l in l1: print(l, file=sys.stderr)  
      for l in l2: print(l, file=sys.stderr)
      sys.exit()
    if currstruc_false or currstruc != stnr:
      print("Invalid SPLIT number:", currstruc, file=sys.stderr)
      print(fnam, file=sys.stderr)
      print("#"+str(stnr), file=sys.stderr)
      for l in l1: print(l, file=sys.stderr)  
      for l in l2: print(l, file=sys.stderr)
      sys.exit()
        
    for lnr, l in enumerate(l1): 
      if lnr == skipline: continue
      print(l)
    for l in l2: print(l)
    
