import sys

import numpy
from math import *

def read_pdb(pdb):
  atoms = []
  lines0 = open(pdb).readlines()
  lines = []
  extralines = []
  for l in lines0:
    if not l.startswith("ATOM"): 
      extralines.append((len(lines), l))
      continue
    x = float(l[30:38])
    y = float(l[38:46])
    z = float(l[46:54])    
    atoms.append((x,y,z))
    lines.append(l)
  return lines, numpy.array(atoms), extralines

def apply_matrix(atoms, pivot, rotmat, trans):
  ret = []  
  for atom in atoms:
    a = atom-pivot
    atom2 = a.dot(rotmat) + pivot + trans
    ret.append(atom2)
  return ret

def fit(atoms1, atoms2):
  # adapted from QKabsch.py by Jason Vertrees. 
  # further adapted from irmsd for fitting by Sjoerd de Vries
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

  # We already have our solution, in the results from SVD.
  # we just need to check for reflections and then produce
  # the rotation.  V and Wt are orthonormal, so their det's
  # are +/-1.0 (and thus products are +/- 1.0 ).
  reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))

  if reflect == -1.0:
	  S[-1] = -S[-1]
	  V[:,-1] = -V[:,-1]
  
  U = V.dot(Wt).transpose()
  RMSD = E0 - (2.0 * sum(S))
  RMSD = numpy.sqrt(abs(RMSD / L))  
  return U, COM1-COM2, RMSD

def cyclesfit(atoms1, atoms2,iterations=5,cutoff=2):
  """
  Performs a PyMol-style fit with rejection cycles
   After every cycle, atoms with more than <cutoff> standard deviations error (in A2)
   get rejected from the selection
  """
  a1, a2 = numpy.array(atoms1), numpy.array(atoms2)
  pivot = numpy.sum(a2,axis=0) / float(len(a2))
  for n in range(iterations):
    rotmat, offset, rmsd = fit(a1,a2)
    aa2 = apply_matrix(a2, pivot, rotmat, offset)
    dif = aa2 - a1
    d = dif.sum(axis=1)
    sd = d * d
    msd = sd.sum()/len(sd)
    #rmsd = sqrt(msd)
    std = numpy.std(sd)
    if (std < 0.1): std = 0.1    
    
    keep = numpy.less(sd,cutoff*std)
    aa1 = a1[keep]
    aa2 = a2[keep]
    
    a1 = aa1
    a2 = aa2    
  return fit(a1, a2)

def select_bb(lines, atoms):
  ret = []
  for l, a in zip(lines, atoms):
    if l[13:15] in ("CA","C ","O ","N "): ret.append(a)
  return ret

def select_ca(lines, atoms):
  ret = []
  for l, a in zip(lines, atoms):
    if l[13:15] in ("CA",): ret.append(a)
  return ret
    
def write_pdb(lines, atoms, extralines):
  count = 0
  pos = 0
  data = zip(lines, atoms)
  while 1:
    while pos < len(extralines):
      p,d = extralines[pos]
      if count < p: break
      print(d.rstrip("\n"))
      pos += 1
    if count == len(data): break
    l,a = data[count]
    ll = l[:30] + "%8.3f%8.3f%8.3f" % (a[0],a[1],a[2]) + l[54:].rstrip("\n")
    print(ll)
    count += 1
  
import sys
import argparse
a = argparse.ArgumentParser(prog="fit.py")
a.add_argument("reference")
a.add_argument("mobile")
a.add_argument("--allatoms", action="store_true")
a.add_argument("--rmsd", action="store_true")
a.add_argument("--iterative", action="store_true")
a.add_argument("--iterative_cycles",type=int,default=5)
a.add_argument("--iterative_cutoff",type=float,default=2)
args = a.parse_args()

#read atoms  
lines1, atoms1, extralines1 = read_pdb(args.reference)
lines2, atoms2, extralines2 = read_pdb(args.mobile)

#select backbone
if args.allatoms:
  atoms1_fit = atoms1
  atoms2_fit = atoms2
else:  
  atoms1_fit = select_bb(lines1, atoms1)
  atoms2_fit = select_bb(lines2, atoms2)
assert len(atoms1_fit) and len(atoms1_fit) == len(atoms2_fit), (len(atoms1_fit), len(atoms2_fit))


if args.iterative:
  #perform a Pymol-style iterative fit
  rotmat, offset, rmsd = cyclesfit(atoms1_fit,atoms2_fit, args.iterative_cycles, args.iterative_cutoff)
else:
  #perform a direct fit
  rotmat, offset, rmsd = fit(atoms1_fit,atoms2_fit)

pivot = numpy.sum(atoms2,axis=0) / float(len(atoms2))
fitted_atoms = apply_matrix(atoms2, pivot, rotmat, offset)
if args.rmsd:
  print "%.3f" % rmsd
else:  
  write_pdb(lines2, fitted_atoms, extralines2)

 
