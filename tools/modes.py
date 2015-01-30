"""
  program  modes
 
  this program calculates the Hinsen modes (Proteins, 1998) using only CA atoms of a protein
  the modes are adapted to the full protein structure by moving each residue as rigid body
  author: Martin Zacharias, Jacobs University Bremen
  Uses the algorithm by Setny et al. to compute normal modes for nucleic acids
  
  Ported to Python by Sjoerd de Vries, 2014-2015
"""

import numpy, sys
from math import *

has_argparse = False
try:
  import argparse  
  has_argparse = True  
except ImportError:
  import optparse  #Python 2.6

#parameters for aacontact
ampl = 1.2
pow2 = 5.0**2

sugar_atoms = ["C1'", "C2'", "C3'", "O4'", "C4'"]


def read_pdb(f):
  coor, res, resn, atom = [], [], [], []
  curr_resid = None
  resindex = 0
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    coor.append((x,y,z))
    resid = l[21:26]
    cresn = l[17:20]
    if resid != curr_resid: 
      curr_resid = resid
      resindex += 1      
    res.append(resindex)
    resn.append(cresn)  
    atom.append(l[12:16])
  return coor, res, resn, atom


if has_argparse:
  parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("pdb",help="PDB file to calculate modes for")
  parser.add_argument("nrmodes",help="number of modes to calculate", type=int, default=20)
else:
  parser = optparse.OptionParser()
  parser.add_argument = parser.add_option
parser.add_argument("--na", "--rna", "--dna",dest="sugar", help="Nucleic acid mode: use the sugar ring in the all-atom PDB file as center instead of the CA atom", action="store_true")
parser.add_argument("--scale", metavar="SCALE CONSTANT", help="Scale all spring force constants", type=float, default=1.0)
parser.add_argument("--aacontact",help="Use the contacts from the all-atom PDB file to determine the spring constants (Setny et al.)", action="store_true")
parser.add_argument("--aapdb",help="all-atom PDB file")

if has_argparse:
  args = parser.parse_args()
else:
  args, positional_args = parser.parse_args()
  args.pdb = None
  args.nrmodes = 20
  if positional_args:
    args.pdb = positional_args[0]
    if len(positional_args) > 1: args.nrmodes = int(positional_args[1])

if args.sugar and not args.aapdb:
  raise ValueError("Nucleic acid mode requires an all-atom PDB reference provided with --aapdb")
if args.aacontact and not args.aapdb:
  raise ValueError("--aacontact requires an all-atom PDB reference provided with --aapdb")

coor, res, resn, atom = read_pdb(args.pdb)
nr_res = len(set(res))

if args.aapdb:
  coor_aa, res_aa, resn_aa, atom_aa = read_pdb(args.aapdb)      
  nr_res_aa = len(set(res_aa))  
  if nr_res != nr_res_aa:
    raise ValueError("All-atom PDB has a different number of residues: %d vs %d" % (nr_res_aa, nr_res))

def make_center(atoms):  
  center = []
  for k in atoms:
    if k.endswith("'") and k[:-1] in atoms:
      continue
    center.append(atoms[k])    
  center = numpy.array(center)  
  center = numpy.mean(center,axis=0)
  return center

centers = []
if args.sugar:
  old_res = None
  for ccoor, cres, cresn, catom in zip(coor_aa, res_aa, resn_aa, atom_aa):
    if cres != old_res:
      if old_res is not None and len(atoms):         
        center = make_center(atoms)
        centers.append(center)                
      atoms = {}
      old_res = cres    
    if catom.strip() in sugar_atoms:       
      atoms[catom.strip()] = ccoor
  if old_res is not None  and len(atoms):
    center = make_center(atoms)
    centers.append(center)        
else:
  for ccoor, catom in zip(coor, atom):
    if catom != " CA ": continue
    centers.append(numpy.array(ccoor))

if len(centers) != nr_res:
  raise ValueError("Cannot determine the residue centers: %d found for %d residues" % (len(centers), nr_res))

if args.aacontact:
  aa = []
  old_res = None
  for ccoor, cres, catom in zip(coor_aa, res_aa, atom_aa):    
    if cres != old_res:
      aa.append([])
      old_res = cres
    if catom.strip().strip("'").startswith("H"): continue      
    aa[-1].append(numpy.array(ccoor))
  assert len(aa) == len(centers)
  
asize = 3*len(centers)  
d2f = numpy.zeros((asize,asize))
springconstants = {}
for i in range(len(centers)):
  for j in range(len(centers)):    
    if i == j: continue
    dx = centers[i] - centers[j]
    r2 = dx.dot(dx)
    if args.aacontact:
      springconst = 0
      for ci in aa[i]:
        for cj in aa[j]:
          dc = ci-cj
          dcsq = dc.dot(dc)
          springconst += ampl * exp(-dcsq/pow2)
    else:      
      springconst = 2 * exp(-r2/16.0)
    springconst *= args.scale
    springconstants[i,j] = springconst
    
    b0 = springconst/r2
    ii = 3 * i
    jj = 3 * j
    
    for k in range(3):
      for l in range(3):
        d2f[ii+k,jj+l] -= b0*dx[k]*dx[l] 
    for k in range(3):
      for l in range(k,3):
        d2f[ii+k,ii+l] += b0*dx[k]*dx[l] 

for i in range(asize):
  for j in range(i): 
    d2f[i,j] = d2f[j,i] 

for i in range(len(centers)):
  for j in range(i):
    springconstants[i,j] = springconstants[j,i]

"""
#print spring constants per coarse-grained atom
for i in range(len(centers)):
  for j in range(len(centers)):
    if i == j: v = 0.0
    else: v = springconstants[i,j]
    print "%11.4f" % v,
  print  
"""

U,s,V = numpy.linalg.svd(d2f)    
ind = numpy.argsort(s)[::-1]
U = U[:, ind]
s = s[ind]
V = V[:, ind]
for i in range(args.nrmodes):  
  n = asize-7-i
  vec = [-float(x) for x in V[n]]
  print " ", i+1, float(s[n])
  for resindex in res:
    pos = 3 * (resindex - 1)
    for n in range(3):
      print vec[pos+n],
    print
    
  print
