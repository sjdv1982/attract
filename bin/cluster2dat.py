import sys
import itertools

def read_next_struc(lines, firstline):

  assert int(firstline[1:]) == 1
  mode = 1
  ret0,ret1 = [],[]
  lenstruc = 1

  for l in lines:
    l = l.rstrip("\n")  
    if mode == 2 and l[0] == "#": 
      mode = 0
      yield ret0,ret1
    if mode == 0:
      assert int(l[1:]) == lenstruc + 1
      mode = 1
      ret0,ret1 = [],[]
      lenstruc += 1
      continue
    if mode == 1:
      if l[:2] == "##":
	ret0.append(l)
      else:
	mode = 2
    if mode == 2:
      ret1.append(l)
  if len(ret0) or len(ret1):
    yield ret0,ret1
  
def read_struc(fil):
  lines = open(fil).xreadlines()
  header = []
  centeredlen = 0
  firstline = None
  for lnr, l in enumerate(lines):
    l = l.rstrip("\n")
    if centeredlen < 2:
      if l.startswith("##"): continue
      header.append(l)
      if l.startswith("#centered"): centeredlen += 1
      continue
    firstline = l
    break
  if firstline is None:
    raise ValueError("Cannot find structures in file %s" % fil)
  return header, read_next_struc(lines,firstline)

cluster = sys.argv[1]
header, structures = read_struc(sys.argv[2])
structures = list(structures)

for h in header: print h

clusterlist = []
for line in open(cluster).readlines():
  tmp = line.split()
  clusterlist.append(int(tmp[3])) #select first representative of cluster
  
stnr = 0
st2nr = 0
for s in structures:
  stnr += 1
  if not stnr in clusterlist:
    continue
  else:
    st2nr += 1
    print "#"+str(st2nr)
    print "##"+str(stnr) + " => cluster"
    l1,l2 = s
    for l in l1: print l
    for l in l2: print l