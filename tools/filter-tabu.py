import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
structures = list(structures)
rmsd_threshold = float(sys.argv[2])
rmsd_files = sys.argv[3:]

def read_rmsd(rmsd_file):
  rmsds = []
  for l in open(rmsd_file).readlines():
    ll = l.split()
    if ll[0] != "RMSD": continue
    rmsds.append([float(v) for v in ll[2:]])
  return rmsds  

all_rmsds = []
for rmsd_file in rmsd_files:
  rmsds = read_rmsd(rmsd_file)
  assert len(rmsds) == len(structures)
  all_rmsds.append(rmsds)

def filter_rmsd(all_rmsds, nr):
  for rmsds in all_rmsds:
    ok = False
    for v in rmsds[nr]:
      if v > rmsd_threshold:
        ok = True
        break
    if not ok: return False
  return True

for h in header: print h
stnr = 0
for r, s in enumerate(structures):
  l1,l2 = s
  if not filter_rmsd(all_rmsds, r): continue
  stnr += 1    
  print "#"+str(stnr)
  print "##"+str(r+1) + " => filter-tabu"
  try:
    for l in l1: print l
    for l in l2: print l
  except IOError:  
    sys.exit()

