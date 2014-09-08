"""
Calculate interface RMSD
usage: python irmsd.py <DAT file> \
 <unbound PDB 1> <bound PDB 1> [<unbound PDB 2> <bound PDB 2>] [...]
 [--allatoms] [--allresidues] [--cutoff <dist cutoff for interface, in A> ]

--allatoms: use all atoms rather than backbone atoms
--allresidues: use also the residues outside the 10 A interface region 

"""
thresh = 10.0
threshsq = thresh * thresh
import sys

import numpy
import collectlibpy as collectlib

def get_interface(boundatoms):
  
  res = []
  pos = 0
  for p in boundatoms:
    cres = []
    ccres = []
    resid = None
    for a,c in zip(p[1],p[0]):
      pos += 1
      cresid = a[21:26]
      #print resid, cresid
      if resid is not None: 
        if cresid != resid:
          cres.append(ccres)
          ccres = []
      ccres.append((pos,a,c))
      resid = cresid   
    cres.append(ccres)
    res.append(cres)
  
  ret = []
  def add_res(l, r, pnr):
    for pos, a, c in r:
      v = (pos, a, c, pnr)
      if v not in l:
        l.append(v)
  def dsq(c1,c2):
    dx = c1[0]-c2[0]
    dy = c1[1]-c2[1]
    dz = c1[2]-c2[2]
    return dx**2+dy**2+dz**2
    
  for pnr,p in enumerate(res):
    ret0 = []
    for r in p:
      ok = False                  

      #check if close to the existing interface
      for pos,a,c1 in r:        
        for pos2,a2,c2,ppnr in ret:
          if ppnr == pnr: continue
          if dsq(c1,c2) < threshsq:
            ok = True
            break
        if ok: break   
      if ok: 
        add_res(ret, r, pnr)
        continue
      
      #if not, check if close to subsequent residues
      for pos,a,c1 in r:        
        for ppnr,pp in enumerate(res[pnr+1:]):
          for rr in pp:     
            for pos2,a2,c2 in rr:
              if dsq(c1,c2) < threshsq:
                ok = True
                break              
            if ok: break     
          if ok: break     
        if ok: break     
      if ok: 
        add_res(ret,r,pnr)
        add_res(ret,rr,pnr+1+ppnr)
    ret += ret0
                  
  ret.sort(key=lambda v: v[0])
  ret = [(v[0],v[1]) for v in ret]
  return ret
  
def get_selection(boundatoms, use_allresidues, use_allatoms):
  
  allatoms = []
  for b in boundatoms: allatoms += b[1]
  if not use_allresidues:
    selatoms = get_interface(boundatoms)
  else:
    selatoms = [(n+1,a) for n,a in enumerate(allatoms)]
      
  if not use_allatoms:
    selatoms = [(n,a) for n,a in selatoms if a[13:15] in ("CA","C ","O ","N ")]
  selected = set([n for n,a in selatoms])

  mask = []
  for n in range(len(allatoms)):
    mask.append((n+1) in selected)

  mask = numpy.array(mask)
  return mask

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


def read_pdb(f):
  ret1 = []
  ret2 = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret1.append((x,y,z))
    ret2.append(l)
  return ret1, ret2
  
if __name__ == "__main__":  
  import os
  ensfiles = []
  modefile = None
  imodefile = None
  name = None
  opt_allatoms = False
  opt_allresidues = False

  anr = 0
  output = None
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
      ensfiles.append((sys.argv[anr+1],sys.argv[anr+2]))
      sys.argv = sys.argv[:anr] + sys.argv[anr+3:]
      anr -= 3
      continue

    if anr <= len(sys.argv)-2 and arg == "--cutoff":
      thresh = sys.argv[anr+1]
			threshsq = thresh * thresh
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 2
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
    
    if anr <= len(sys.argv)-2 and arg == "--thresh":
      thresh = float(sys.argv[anr+1])
      threshsq = thresh*thresh
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

  unbounds = []
  bounds = []

  
  for n in range(2, len(sys.argv), 2):
    unbounds.append(sys.argv[n])
    bounds.append(sys.argv[n+1])

  if len(bounds) == 1 and opt_allresidues == False:
    raise Exception("Cannot determine the interface for a single PDB")

  initargs = [sys.argv[1]] + unbounds
  if modefile: initargs += ["--modes", modefile]
  if imodefile: initargs += ["--imodes", imodefile]
  for nr, ensfile in ensfiles:
    initargs += ["--ens", nr, ensfile]

  collectlib.collect_init(initargs)

  boundatoms = []
  for b in bounds:
    boundatoms.append(read_pdb(b))

  boundsizes = [len(b[1]) for b in boundatoms]
  unboundsizes = []
  start = 0
  for i in collectlib.ieins[:len(unbounds)]:
    unboundsizes.append(i-start)
    start = i

  for bname, ubname, bsize, ubsize in zip(bounds,unbounds,boundsizes,unboundsizes):
    if bsize != ubsize:
      raise Exception("Different atom numbers: %s: %d, %s: %d" % (ubname, ubsize, bname, bsize))

  allboundatoms = []
  for b in boundatoms: allboundatoms += b[0]
  sel = get_selection(boundatoms, opt_allresidues, opt_allatoms)
  allboundatoms = numpy.array(allboundatoms)
  fboundatoms = allboundatoms[sel]

  nstruc = 0
  f1 = sys.stdout
  if output is not None:
    f1 = open(output,'w')
  while 1:
    sys.stdout.flush()
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
    fcoor = coor[sel]
    f1.write(str(nstruc)+" %.3f\n" % irmsd(fboundatoms,fcoor))
    #for co in coor: print co[0], co[-1]
