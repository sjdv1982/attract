import sys
rpdb = sys.argv[1]
atomsfile = sys.argv[2]

pdb = []
for l in open(rpdb).readlines():
  if l.startswith("ATOM"): pdb.append(l)
  
count = 0  
for l in open(atomsfile).readlines():
  count += 1
  print "MODEL", count
  ll = l.split()
  n1 = int(ll[0])
  n2 = int(ll[1])
  a1 = pdb[n1-1]
  a2 = pdb[n2-1]
  a1 = a1[:21] + chr(64+count) + a1[22:]
  a2 = a2[:21] + chr(64+count) + a2[22:]
  print a1,
  print a2,  
  print "ENDMDL"
  if count == 26: count = 0
