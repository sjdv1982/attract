"""
Aligns protein complexes on the interface region of the bound PDB
Does not align the monomers (use align-monomers-interface.py for that)
Writes out a DAT file 
usage: python align-interface.py <DAT file> \
 <unbound PDB 1> <bound PDB 1> [<unbound PDB 2> <bound PDB 2>] [...]
 [--allatoms] [--allresidues]

--allatoms: use all atoms rather than backbone atoms
--allresidues: use also the residues outside the 10 A interface region 

"""
thresh = 10.0
threshsq = thresh * thresh
import sys

import numpy
import collectlibpy as collectlib
from _read_struc import read_struc
from euler2rotmat import euler2rotmat
from rotmat2euler import rotmat2euler

def align_interface(atoms1, atoms2):
  # adapted from QKabsch.py by Jason Vertrees. 
  L = len(atoms1)
  assert( L > 0 )

  # must alway center the two proteins to avoid
  # affine transformations.  Center the two proteins
  # to their selections.
  COM = numpy.sum(atoms2,axis=0) / float(L)
  atoms2 = atoms2 - COM

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

  if reflect < -0.99:
	  S[-1] = -S[-1]
	  V[:,-1] = -V[:,-1]
  rotmatd = numpy.dot( V, Wt)
  rotmatd = numpy.transpose(rotmatd)
  rotmatd = [[min(xx,1) for xx in x] for x in rotmatd]
  rotmatd = numpy.array([[max(xx,-1) for xx in x] for x in rotmatd])
  return rotmatd, COM

import sys, os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib

if __name__ == "__main__":  
  ensfiles = []
  ensembles = []
  modefile = None
  imodefile = None
  opt_allatoms = False
  opt_allresidues = False

  anr = 0
  output = None
  atomnames = ("CA","C","O","N")
  while 1:
    anr += 1
        
    if anr > len(sys.argv)-1: break  
    arg = sys.argv[anr]

    if arg == "--allatoms": 
      sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
      opt_allatoms = True
      anr -= 1
      continue
    
    if arg == "--allresidues": 
      sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
      opt_allresidues = True
      anr -= 1
      continue  
    
    
    if anr <= len(sys.argv)-3 and arg == "--ens":
      ens = sys.argv[anr+1]
      ensfiles.append((ens,sys.argv[anr+2]))
      ensembles.append(ens)
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
    
    if anr <= len(sys.argv)-2 and arg == "--output":
      output = sys.argv[anr+1]
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
      continue
    if arg.startswith("--"): raise Exception("Unknown option '%s'" % arg)
      

  if len(sys.argv) < 4 or len(sys.argv) % 2:
    raise Exception("Please supply an even number of PDB files (unbound, bound)")
  
  unboundfiles = []
  boundfiles = []
  for n in range(2, len(sys.argv), 2):
    unboundfiles.append(sys.argv[n])
    boundfiles.append(sys.argv[n+1])

  if len(boundfiles) == 1 and opt_allresidues == False:
    raise Exception("Cannot determine the interface for a single PDB")

  bounds = [rmsdlib.read_pdb(f) for f in boundfiles]
  unbounds = [rmsdlib.read_pdb(f) for f in unboundfiles]

  struc_header,structures = read_struc(sys.argv[1])
  pivots = []
  for hnr,h in enumerate(struc_header):
    if not h.startswith("#pivot"): continue
    hh = h.split()
    assert len(hh) == 5 and hh[1] == str(hnr+1), h
    pivot = numpy.array([float(v) for v in hh[2:5]])
    pivots.append(pivot)
  
  initargs = [sys.argv[1]] + unboundfiles
  if modefile: initargs += ["--modes", modefile]
  if imodefile: initargs += ["--imodes", imodefile]
  for nr, ensfile in ensfiles:
    initargs += ["--ens", nr, ensfile]

  collectlib.collect_init(initargs)
  unboundsizes = [len(list(p.atoms())) for p in unbounds]
  collectlib.check_sizes(unboundsizes, unboundfiles)

  rmsdlib.check_pdbs(unbounds, bounds)

  allboundatoms = []
  for p in bounds:
    for c in p.coordinates():
      allboundatoms.append(c)
  allboundatoms = numpy.array(allboundatoms)

  sel = numpy.array([True] * len(allboundatoms))
  if not opt_allresidues:
    imask = rmsdlib.build_interfacemask(bounds, thresh)
    sel = (sel & imask)
  if not opt_allatoms:
    amask = rmsdlib.build_atommask(bounds, atomnames)
    sel = (sel & amask)
  allboundatoms = numpy.array(allboundatoms)
  fboundatoms = numpy.array(allboundatoms[sel])
  icom = numpy.sum(fboundatoms,axis=0) / float(len(fboundatoms))
  fboundatoms = fboundatoms - icom
  
  nstruc = 0
  f1 = sys.stdout
  if output is not None:
    f1 = open(output,'w')
  print >> f1, "\n".join(struc_header)  
  for structure in structures:
    ligands = structure[1]
    assert len(ligands) == len(pivots), (len(ligands),  len(pivots))
    sys.stdout.flush()
    result = collectlib.collect_next()
    if result: break
    nstruc += 1
    coor = collectlib.collect_all_coor()
    coor = numpy.array(coor)
    fcoor = coor[sel]
    irotmat, ipivot = align_interface(fboundatoms,fcoor) 
    
    print >> f1, "#" + str(nstruc)
    if len(structure[0]):
      print >> f1, "\n".join(structure[0])
    for lignr, lig in enumerate(ligands):
      ll = lig.split()      
      if str(lignr+1) in ensembles:
        ens = ll[0] + " "
        dofs = ll[1:]
      else:
        ens = ""
        dofs = ll
      rotdofs = [float(d) for d in dofs[:3]]        
      transdofs = numpy.array([float(d) for d in dofs[3:6]])
      rotmat = numpy.array(euler2rotmat(*rotdofs))
      newrotmat = numpy.dot(irotmat.T, rotmat)      
      
      #rotate around ipivot
      transdofs_global = transdofs + pivots[lignr] #for rotation around the global origin
      transdofs_centered = transdofs_global - ipivot #for rotation around ipivot
      transdofs_new = numpy.dot(irotmat.T, transdofs_centered) - pivots[lignr]
      
      newtransdofs = transdofs_new + icom
      newrotdofs = rotmat2euler(newrotmat)  
      lignew = "  " + ens + "%.6f %.6f %.6f " % tuple(newrotdofs) + "%.6f %.6f %.6f " % tuple(newtransdofs) + " ".join(dofs[6:])
      print >> f1, lignew    
    
