import sys
chain = 'A'
for l in open(sys.argv[1]):
  if not l.startswith("ATOM"):
    if l.startswith("TER"):
      chain = chr(ord(chain)+1)
    print l
    continue
  print l[:21] + chain + l[22:],