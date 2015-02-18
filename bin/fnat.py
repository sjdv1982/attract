import sys
import numpy
import collectlibpy as collectlib
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib
from scipy.spatial.distance import cdist

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

unboundfiles = []
boundfiles = []
for n in range(3, len(sys.argv), 2):
  unboundfiles.append(sys.argv[n])
  boundfiles.append(sys.argv[n+1])

if len(boundfiles) == 1 :
  raise Exception("Cannot determine the contacts for a single PDB")

bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
unbounds = [rmsdlib.read_pdb(f) for f in unboundfiles]

initargs = [sys.argv[1]] + unboundfiles
if modefile: initargs += ["--modes", modefile]
if imodefile: initargs += ["--imodes", imodefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)
unboundsizes = [len(list(p.atoms())) for p in unbounds]
collectlib.check_sizes(unboundsizes, unboundfiles)

unbound_hmask = rmsdlib.build_hydrogenmask(unbounds)
for p in unbounds + bounds: 
  p.remove_hydrogens()

rmsdlib.check_pdbs(unbounds, bounds)

contacts = rmsdlib.build_contactfilters(bounds, cutoff)

nstruc = 0
f1 = sys.stdout
if output is not None:
  f1 = open(output,'w')
  
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
  coor = numpy.compress(unbound_hmask, coor, axis=0)
  fcount = 0
  for filter1, filter2 in contacts:
    coor1 = numpy.take(coor, filter1, axis=0)
    coor2 = numpy.take(coor, filter2, axis=0)
    distances = cdist(coor1, coor2,'sqeuclidean')
    mindistance = distances.min()
    if mindistance < cutoffsq:
      fcount += 1
  f1.write("%.2f" % (float(fcount)/len(contacts)) +'\n')
    
