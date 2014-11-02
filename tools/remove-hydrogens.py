import sys
for l in open(sys.argv[1]):
  if l.startswith("ATOM") or l.startswith("HETATM"):
    code = l[12:16].strip()
    if code.startswith("H"): 
      continue
  print l,