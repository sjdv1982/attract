#Fits two PDB files, describes the euler angles and translation 

import sys
from math import *
from rotmat2euler import rotmat2euler

import numpy

def euler(atoms1, atoms2):
  # adapted from irmsd by Sjoerd de Vries. 
  # adapted from QKabsch.py by Jason Vertrees. 
  L = len(atoms1)
  assert( L > 0 )

  # must alway center the two proteins to avoid
  # affine transformations.  Center the two proteins
  # to their selections.
  COM1 = numpy.sum(atoms1,axis=0) / float(L)
  COM2 = numpy.sum(atoms2,axis=0) / float(L)
  dx = COM2[0]-COM1[0]
  dy = COM2[1]-COM1[1]
  dz = COM2[2]-COM1[2]
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
  #print V
  #print Wt
  
  rotmatd = numpy.dot( V, Wt)
  rotmatd = numpy.transpose(rotmatd)
  rotmatd = [[min(xx,1) for xx in x] for x in rotmatd]
  rotmatd = [[max(xx,-1) for xx in x] for x in rotmatd]
  #print ["%.3f" % v for v in rotmatd[0]]
  #print ["%.3f" % v for v in rotmatd[1]]
  #print ["%.3f" % v for v in rotmatd[2]]

  phi,ssi,rot = rotmat2euler(rotmatd)
  return phi, ssi, rot, dx, dy, dz


def read_pdb(f):
  ret = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    if l[13] == 'H': continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret.append((x,y,z))
  return ret

if __name__ == "__main__":  

  if len(sys.argv) !=  3:
    raise Exception("Please supply two PDB files")

  atoms1 = read_pdb(sys.argv[1])
  atoms2 = read_pdb(sys.argv[2])

  if len(atoms1) != len(atoms2):
    raise Exception("Different atom numbers: %s: %d, %s: %d" % (sys.argv[1], len(atoms1), sys.argv[2], len(atoms2)))

  atoms1 = numpy.array(atoms1)
  atoms2 = numpy.array(atoms2)

  ret = euler(atoms1, atoms2)
  ret = ["%.5f" % x if fabs(x) > 1e-5 else 0 for x in ret ]
  phi, ssi, rot, dx, dy, dz = ret
  print phi, ssi, rot, dx, dy, dz
  