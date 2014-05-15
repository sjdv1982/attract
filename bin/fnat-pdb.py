import sys

import numpy

ensfiles = []
modefile = None
imodefile = None
output = None
anr = 0
if len(sys.argv) < 5 or not len(sys.argv) % 2:
  raise Exception("Please supply a multi-model PDB file, a cutoff and a number of bound PDB files")

cutoff = float(sys.argv[2])
cutoffsq = cutoff*cutoff

bounds = []

def read_multi_pdb(f):  
  endmodel = False
  allcoor = [[]]
  coor = allcoor[-1]
  for l in open(f):
    if l.startswith("MODEL"):      
      allcoor = [[]]
      coor = allcoor[-1]
    if l.startswith("TER"):
      allcoor.append([])
      coor = allcoor[-1]
    if l.startswith("ENDMDL"):      
      endmodel = True
      allcoor = [coor for coor in allcoor if len(coor)]
      yield allcoor     
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    coor.append((x,y,z))
    endmodel = False
  allcoor = [coor for coor in allcoor if len(coor)]
  if not endmodel: yield allcoor

def read_pdb(f):
  ret, res, resnam = [], [],{}
  curr_resid = None
  resindex = 0
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret.append((x,y,z))
    resid = l[21:26]
    if resid != curr_resid: 
      curr_resid = resid
      resindex += 1
      resnam[resindex] = resid
    res.append(resindex)  
  return ret, res, resnam
  
bounds = sys.argv[3:]

boundatoms = []
boundres = []
boundresnam = []
boundsizes = []
for b in bounds:
  atoms, res, resnam = read_pdb(b)
  boundatoms.append(atoms)
  boundsizes.append(len(atoms))
  boundres.append(res)
  boundresnam.append(resnam)

  
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
  
start = 0
for i in read_multi_pdb(sys.argv[1]):
  if len(i) != len(boundsizes):
      raise Exception(
"Chain size difference between PDB and templates: %d vs %d chains, %s vs %s residues" % (len(i), len(boundsizes), [len(ii) for ii in i], boundsizes)
    )
  for chainnr, chain in enumerate(i):
    if len(chain) != boundsizes[chainnr]:
      raise Exception(
"Atom size difference between PDB and template: PDB %s: %d vs %d atoms" % (bounds[chainnr], len(chain), boundsizes[chainnr])
    )
  coor = []
  for chain in i: coor += chain
  coor = numpy.array(coor)
  nstruc += 1
  fcount = 0
  for maskA, maskB in contacts:
    coorA = coor[maskA]
    coorB = coor[maskB]
    distances = cdist(coorA, coorB,'sqeuclidean')
    mindistance = distances.min()
    #print(mindistance)
    if mindistance < cutoffsq:
      fcount += 1
  f.write(str(float(fcount)/len(contacts))+'\n')
    
