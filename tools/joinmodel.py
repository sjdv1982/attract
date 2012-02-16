import sys

for fnr,f in enumerate(sys.argv[1:]):
  print("MODEL %d" % (fnr + 1))
  for l in open(f):
    if l.startswith("END"): continue
    print(l.rstrip("\n"))
  print("ENDMDL")
