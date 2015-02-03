import sys

import numpy
import collectlibpy as collectlib
import os
ensfiles = []
modefile = None
imodefile = None
output = None
name = None
anr = 0
while 1:
  anr += 1

  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]
  
  if anr <= len(sys.argv)-3 and arg == "--ens":
    ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
    sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
    anr -= 3
    continue

  if anr <= len(sys.argv)-2 and arg == "--modes":
    modefile = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
  
  if anr <= len(sys.argv)-2 and arg == "--imodes":
    imodefile = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
# Suppport for direct IRMSD with iATTRACT files
  if anr <= len(sys.argv)-2 and arg == "--name":
    name = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
    
  if anr <= len(sys.argv)-2 and arg == "--output":
    output = sys.argv[anr+1]
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 2
    continue
if len(sys.argv) < 5 or not len(sys.argv) % 2:
  raise Exception("Please supply a .DAT file, a cutoff and an even number of PDB files (unbound, bound)")

cutoff = float(sys.argv[2])
cutoffsq = cutoff*cutoff

unbounds = []
bounds = []

def read_pdb(f):
  ret, res, resnam = [], [],{}
  curr_resid = None
  resindex = 0
  for l in open(f):
    if not l.startswith("ATOM"): continue
    atomcode = l[12:16].strip()
    if atomcode.startswith("H"): continue  
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret.append((x,y,z))
    resid = l[21:26]
    if resid != curr_resid: 
      curr_resid = resid
      resindex += 1
      resnam[resindex] = resid
    res.append(resindex)  
  return ret, res, resnam
  
for n in range(3, len(sys.argv), 2):
  unbounds.append(sys.argv[n])
  bounds.append(sys.argv[n+1])

initargs = [sys.argv[1]] + unbounds
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)

boundatoms = []
boundres = []
boundresnam = []
for b in bounds:
  atoms, res, resnam = read_pdb(b)
  boundatoms.append(atoms)
  boundres.append(res)
  boundresnam.append(resnam)

unboundatoms = []
unboundres = []
unboundresnam = []
for ub in unbounds:
  atoms, res, resnam = read_pdb(ub)
  unboundatoms.append(atoms)
  unboundres.append(res)
  unboundresnam.append(resnam)
  
boundsizes = [len(b) for b in boundatoms]
unboundsizes = [len(ub) for ub in unboundatoms]
start = 0
for inr,i in enumerate(collectlib.ieins[:len(unbounds)]):
  collectsize = i-start
  if collectsize != unboundsizes[inr]:
    raise Exception(
"Parsing difference between collect and Python: PDB %s: %d vs %d atoms" % (unbounds[inr], collectsize, unboundsizes[inr])
)
  start = i

for bname, ubname, bres, ubres, bresnam, ubresnam in \
 zip(bounds,unbounds,boundres,unboundres, boundresnam, unboundresnam):
  if len(bres) != len(ubres):
    print("ERROR: Different number of atoms: %s: %d, %s: %d" % (ubname, len(ubres), bname, len(bres)))
  if bres != ubres:
    for pos,item in enumerate(zip(bres, ubres)):
      rb, rub = item
      rbnam = bresnam[rb]
      rubnam = bresnam[rub]
      if rb != rub:
        raise Exception("Residue layout differs: %s and %s, atom %d: '%s'(residue %d) - '%s'(residue %d) " 
         % (ubname, bname, pos+1, rbnam, rb, rubnam, rub))


#print("START")
p = []
masks = []
for partner in range(len(bounds)):
  coor = numpy.array(boundatoms[partner])
  res = numpy.array(boundres[partner])
  resnr = res[-1]
  p.append((coor, res, resnr))
  pmasks = []
  for resid in range(resnr):  
    mask = numpy.ma.masked_equal(res,resid+1).mask
    coorres = coor[mask]
    pmasks.append((mask, coorres))
  masks.append(pmasks)

#build bound-form contacts
contacts0 = []
from scipy.spatial.distance import cdist
for partner1 in range(len(bounds)):
  coor1, res1, resnr1 = p[partner1]
  for partner2 in range(partner1+1, len(bounds)):
    coor2, res2, resnr2 = p[partner2]
    for resid1 in range(resnr1):
      mask1, coor1res = masks[partner1][resid1]
      for resid2 in range(resnr2):
        mask2, coor2res = masks[partner2][resid2]
        distances = cdist(coor1res, coor2res,'sqeuclidean')
        mindistance = distances.min()
        if mindistance < cutoffsq:
          contacts0.append((partner1, mask1, partner2, mask2))

#build contacts as selection masks
contacts = []
allatoms = 0
offsets = []
for atoms in boundatoms:
  offsets.append(allatoms)
  allatoms += len(atoms)
offsets.append(allatoms)  
mask0 = numpy.array([False] * allatoms)
for partner1, mask1, partner2, mask2 in contacts0:
  maskA = mask0.copy()
  maskB = mask0.copy()
  maskA[offsets[partner1]:offsets[partner1+1]] = mask1
  maskB[offsets[partner2]:offsets[partner2+1]] = mask2
  contacts.append((maskA, maskB))
del contacts0
  
nstruc = 0
f = sys.stdout
if output is not None:
  f = open(output,'w')
  
while 1:
  if name is not None: 
    newargs = initargs + ['--imodes','flexm-'+str(nstruc+1)+name+'.dat']
    if not os.path.exists('flexm-'+str(nstruc+1)+name+'.dat'):
      break
    collectlib.collect_iattract(newargs)
    
  result = collectlib.collect_next()
  if result: break
  nstruc += 1
  coor = collectlib.collect_all_coor()
  coor = numpy.array(coor)
  fcount = 0
  for maskA, maskB in contacts:
    coorA = coor[maskA]
    coorB = coor[maskB]
    distances = cdist(coorA, coorB,'sqeuclidean')
    mindistance = distances.min()
    #print(mindistance)
    if mindistance < cutoffsq:
      fcount += 1
  f.write("%.2f" % (float(fcount)/len(contacts)) +'\n')
    
