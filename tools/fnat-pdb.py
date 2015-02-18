import sys
import numpy
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib
from scipy.spatial.distance import cdist

ensfiles = []
modefile = None
imodefile = None
output = None
anr = 0
if len(sys.argv) < 4:
  raise Exception("Please supply a multi-model PDB file, a cutoff and a number of bound PDB files")

cutoff = float(sys.argv[2])
cutoffsq = cutoff*cutoff

boundfiles = sys.argv[3:]

if len(boundfiles) == 1 :
  raise Exception("Cannot determine the contacts for a single PDB")

bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
for p in bounds: 
  p.remove_hydrogens()

contacts = rmsdlib.build_contactfilters(bounds, cutoff)

nstruc = 0
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')

for unbounds in rmsdlib.read_multi_pdb(sys.argv[1]):
  coor = []  
  for p in unbounds: 
    p.remove_hydrogens()
    for c in p.coordinates():
      coor.append(c)
  coor = numpy.array(coor)
  rmsdlib.check_pdbs(unbounds, bounds)
  nstruc += 1
  fcount = 0
  for filter1, filter2 in contacts:
    coor1 = numpy.take(coor, filter1, axis=0)
    coor2 = numpy.take(coor, filter2, axis=0)
    distances = cdist(coor1, coor2,'sqeuclidean')
    mindistance = distances.min()
    if mindistance < cutoffsq:
      fcount += 1
  f1.write("%.2f" % (float(fcount)/len(contacts)) +'\n')
