import sys

import numpy
import collectlibpy as collectlib

def irmsd(atoms1, atoms2):
  # adapted from QKabsch.py by Jason Vertrees. 
  L = len(atoms1)
  assert( L > 0 )

  # must alway center the two proteins to avoid
  # affine transformations.  Center the two proteins
  # to their selections.
  COM1 = numpy.sum(atoms1,axis=0) / float(L)
  COM2 = numpy.sum(atoms2,axis=0) / float(L)
  atoms1 = atoms1 - COM1
  atoms2 = atoms2 - COM2

  # Initial residual, see Kabsch.
  E0 = numpy.sum( numpy.sum(atoms1 * atoms1,axis=0),axis=0) + numpy.sum( numpy.sum(atoms2 * atoms2,axis=0),axis=0)

  #
  # This beautiful step provides the answer.  V and Wt are the orthonormal
  # bases that when multiplied by each other give us the rotation matrix, U.
  # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
  V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(atoms1), atoms2))

  # HACK COMMENT:
  # Numpy has some strangeness with returning floats from
  # its calculation.  So, I made a rather silly work around,
  # but, it seems to work!  I cast the float to a string,
  # then back to a float; this works.  See below.

  # We already have our solution, in the results from SVD.
  # we just need to check for reflections and then produce
  # the rotation.  V and Wt are orthonormal, so their det's
  # are +/-1.0 (and thus products are +/- 1.0 ).
  reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))

  if reflect == -1.0:
	  S[-1] = -S[-1]
	  V[:,-1] = -V[:,-1]

  RMSD = E0 - (2.0 * sum(S))
  RMSD = numpy.sqrt(abs(RMSD / L))
  return RMSD

ensfiles = []
modefile = None
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

if len(sys.argv) < 4 or len(sys.argv) % 2:
  raise Exception("Please supply an even number of PDB files (unbound, bound)")

unbounds = []
bounds = []

def read_pdb(f):
  ret = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret.append((x,y,z))
  return ret
  
for n in range(2, len(sys.argv), 2):
  unbounds.append(sys.argv[n])
  bounds.append(sys.argv[n+1])

initargs = [sys.argv[1]] + unbounds
if modefile: initargs += ["--modes", modefile]
for nr, ensfile in ensfiles:
  initargs += ["--ens", nr, ensfile]

collectlib.collect_init(initargs)

boundatoms = []
for b in bounds:
  boundatoms.append(read_pdb(b))

boundsizes = [len(b) for b in boundatoms]
unboundsizes = []
start = 0
for i in collectlib.ieins[:len(unbounds)]:
  unboundsizes.append(i-start)
  start = i

for bname, ubname, bsize, ubsize in zip(bounds,unbounds,boundsizes,unboundsizes):
  if bsize != ubsize:
    raise Exception("Different atom numbers: %s: %d, %s: %d" % (ubname, ubsize, bname, bsize))

allboundatoms = []
for b in boundatoms: allboundatoms += b
allboundatoms = numpy.array(allboundatoms)

nstruc = 0
while 1:
  result = collectlib.collect_next()
  if result: break
  nstruc += 1
  coor = collectlib.collect_all_coor()
  coor = numpy.array(coor)
  print nstruc, "%.3f" % irmsd(allboundatoms,coor)
  #for co in coor: print co[0], co[-1]
