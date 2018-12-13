import os
import sys

class Atom(object):
  def __init__(self , l):
    self.name = l[12:16].strip()
    self.x = float(l[30:38])
    self.y = float(l[38:46])
    self.z = float(l[46:54])
    self.resid = l[21:27]
    self.chain = l[21] #readonly
    self.resnr = int(l[22:26]) #readonly
    self.resname = l[17:20].strip()
    self.line = l
  def write(self, filehandle):
    l = self.line
    o = l[:12]
    if len(self.name) == 4:
      o += "%4s " % self.name
    else:
      o += " %-3.3s " % self.name
    o += "%-3.3s" % self.resname
    o += l[20]
    o += self.resid
    o += l[27:30]
    o += "%8.3f" % self.x
    o += "%8.3f" % self.y
    o += "%8.3f" % self.z
    o += l[54:]
    filehandle.write(o)

class PDB(object):
  def __init__(self, filename):
    self.filename = filename
    self._residues = {}
    self._res = []
    self._lastres = None
    self._lastresid = None
  def _add_atom(self, a):
    i = a.resid
    i2 = self._lastresid
    if i != self._lastres:
      if i not in self._residues:
        self._res.append(i)
        self._residues[i] = []
        i2 = i
      else: #later residue with the same resid
        count = 2
        while 1:
          i2 = (i, count)
          if i2 not in self._residues: break
          count += 1
        self._res.append(i2)
        self._residues[i2] = []
    self._residues[i2].append(a)
    self._lastres = i
    self._lastresid = i2
  def residues(self):
    for n in self._res:
      yield self._residues[n]
  def atoms(self):
    for n in self._res:
      for a in self._residues[n]:
        yield a
  def coordinates(self):
    coors = []
    for n in self._res:
      for a in self._residues[n]:
        coors.append((a.x,a.y,a.z))
    return coors

  def remove_hydrogens(self):
    for n in list(self._res):
      res = []
      for a in self._residues[n]:
        if a.name[0] != "H": res.append(a)
      if not res:
        self._residues.pop(n)
        self._res.remove(n)
      else:
        self._residues[n] = res

  def select_atoms(self, atomnames):
    for n in list(self._res):
      res = []
      for a in self._residues[n]:
        if a.name in atomnames: res.append(a)
      if not res:
        self._residues.pop(n)
        self._res.remove(n)
      else:
        self._residues[n] = res

  def write(self, fil):
    filehandle = open(fil, "w")
    for n in self._res:
      for a in self._residues[n]:
        a.write(filehandle)
    filehandle.close()

def read_pdb(f):
  assert os.path.exists(f), f
  p = PDB(f)
  for l in open(f):
    if not l.startswith("ATOM"): continue
    a = Atom(l)
    p._add_atom(a)
  return p

def read_multi_pdb(f):
  count = 1
  chain = 1
  endmodel = False
  p = PDB("model %d chain %d" % (count, chain))
  pdbs = [p]
  for l in open(f):
    if l.startswith("MODEL"):
      chain = 1
      if endmodel:
        count += 1
        p = PDB("model %d chain %d" % (count, chain))
        pdbs = [p]
    elif l.startswith("TER"):
      chain += 1
      p = PDB("model %d chain %d" % (count, chain))
      pdbs.append(p)
    elif l.startswith("ENDMDL") and not endmodel:
      endmodel = True
      if not len(pdbs[-1]._res): pdbs = pdbs[:-1]
      yield pdbs
    if not l.startswith("ATOM"): continue
    endmodel = False
    a = Atom(l)
    p._add_atom(a)
  if not endmodel:
    if not len(pdbs[-1]._res): pdbs = pdbs[:-1]
    yield pdbs

def read_multi_attract_pdb(f):
  partner = 1
  p = PDB("partner %d" % (partner))
  pdbs = [p]
  partners = None
  for l in open(f):
    if partners is None:
      partners = int(l)
      continue
    if l.startswith("TER"):
      if not len(p._res):
        pdbs = pdbs[:-1]
      partner += 1
      p = PDB("partner %d" % (partner))
      pdbs.append(p)
    if not l.startswith("ATOM"): continue
    a = Atom(l)
    p._add_atom(a)
  if not len(p._res):
    pdbs = pdbs[:-1]
  assert len(pdbs) == partners, (len(pdbs), partners)
  return pdbs

def check_pdbs(unbounds, bounds):
  for ub, b in zip(unbounds, bounds):
    ubha = [a for a in ub.atoms()]
    bha = [a for a in b.atoms()]
    if len(ubha) != len(bha):
      raise Exception("Different number of atoms: %s: %d, %s: %d" % (ub.filename, len(ubha), b.filename, len(bha)))
    atomcount = 0
    for rub0, rb0 in zip(ub.residues(), b.residues()):
      rub = [a for a in rub0]
      rb = [a for a in rb0]
      if len(rub) != len(rb):
        raise Exception("Residue layout differs: %s, residue %s: %d atoms; %s, residue %s: %d atoms" % \
         (ub.filename, rub[0].resid.strip(), len(rub), b.filename, rb[0].resid.strip(), len(rb)))
      for aub, ab in zip(rub, rb):
        atomcount += 1
        if aub.name != ab.name or aub.resname != ab.resname:
          raise Exception("Residue layout differs, heavy atom %d: %s, residue %s: %s; %s, residue %s: %s" % \
           (atomcount+1, ub.filename, rub[0].resid.strip(), aub.resname+" "+aub.name, b.filename, rb[0].resid.strip(), ab.resname+" "+ab.name))

def build_hydrogenmask(pdbs):
  import numpy
  hmask = []
  for pdb in pdbs:
    for a in pdb.atoms():
      hmask.append(a.name[0] != "H")
  return numpy.array(hmask)

def build_atommask(pdbs, atoms):
  import numpy
  amask = []
  for pdbnr, pdb in enumerate(pdbs):
    ok = False
    for a in pdb.atoms():
      found = a.name in atoms
      amask.append(found)
      if found:
          ok = True
    if not ok:
        raise Exception("Atommask: zero atoms selected for molecule %d" % (pdbnr+1))
  return numpy.array(amask)

def _build_resfilters(pdbs):
  import numpy
  pp = []
  #build coordinates and residue filters
  offset = 0
  for p in pdbs:
    coor = numpy.array(p.coordinates())
    resfilters = {}
    for r in p._residues:
      resfilter = [anr for anr,a in enumerate(p.atoms()) if a.resid == r]
      resfilters[r] = numpy.array(resfilter)
    pp.append((coor, resfilters, offset))
    offset += len(coor)
  return pp

def build_interfacemask(pdbs, cutoff):
  import numpy
  cutoff=float(cutoff)
  cutoffsq = cutoff*cutoff

  pp = _build_resfilters(pdbs)

  totlen = sum([len(list(p.atoms())) for p in pdbs])
  imask = [False] * totlen

  from scipy.spatial.distance import cdist
  for partner1 in range(len(pdbs)):
    coor1, resfilters, offset = pp[partner1]
    for resfilter in resfilters.values():
      coor1res = numpy.take(coor1,resfilter,axis=0)
      for partner2 in range(len(pdbs)):
        if partner1 == partner2: continue
        coor2 = pp[partner2][0]
        distances = cdist(coor1res, coor2,'sqeuclidean')
        mindistance = distances.min()
        if mindistance < cutoffsq:
          filter = resfilter+offset
          for a in filter:
            imask[a] = True
          break
  return numpy.array(imask)

def build_contactfilters(pdbs, cutoff):
  import numpy
  cutoffsq = cutoff*cutoff

  pp = _build_resfilters(pdbs)

  ret = []
  from scipy.spatial.distance import cdist
  for partner1 in range(len(pdbs)):
    coor1, resfilters1, offset1 = pp[partner1]
    for resfilter1 in resfilters1.values():
      coor1res = numpy.take(coor1,resfilter1,axis=0)
      filter1 = resfilter1+offset1
      for partner2 in range(partner1+1, len(pdbs)):
        coor2, resfilters2, offset2 = pp[partner2]
        for resfilter2 in resfilters2.values():
          coor2res = numpy.take(coor2,resfilter2,axis=0)
          distances = cdist(coor1res, coor2res,'sqeuclidean')
          mindistance = distances.min()
          if mindistance < cutoffsq:
            filter2 = resfilter2+offset2
            ret.append((filter1, filter2))
  return ret

def fit(atoms1, atoms2):
  import numpy
  # adapted from QKabsch.py by Jason Vertrees.
  # further adapted from irmsd for fitting by Sjoerd de Vries
  assert len(atoms1) == len(atoms2)
  assert len(atoms1) > 0
  L = len(atoms1)

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

  if reflect < -0.999:
          S[-1] = -S[-1]
          V[:,-1] = -V[:,-1]

  U = V.dot(Wt).transpose()
  RMSD = E0 - (2.0 * sum(S))
  RMSD = numpy.sqrt(abs(RMSD / L))
  return U, COM1-COM2, RMSD


def multifit(array_atoms1, atoms2):
  """
  Fits an array of atom sets (array_atoms1) onto an atom set (atoms2)
  """
  import numpy
  assert isinstance(array_atoms1, numpy.ndarray)
  assert isinstance(atoms2, numpy.ndarray)

  assert len(array_atoms1.shape) == 3
  assert len(atoms2.shape) == 2

  assert len(atoms2) > 0
  assert array_atoms1.shape[2] == atoms2.shape[1] == 3
  assert len(atoms2) == array_atoms1.shape[1], (atoms2.shape, array_atoms1.shape)
  L = len(atoms2)

  # must alway center the two proteins to avoid
  # affine transformations.  Center the two proteins
  # to their selections.
  COM1 = numpy.sum(array_atoms1,axis=1) / float(L)
  COM2 = numpy.sum(atoms2,axis=0) / float(L)

  array_atoms1 = array_atoms1 - COM1[:,numpy.newaxis, :]
  atoms2 = atoms2 - COM2

  # Initial residual, see Kabsch.
  E0 = numpy.sum( numpy.sum(array_atoms1 * array_atoms1,axis=1),axis=1) + numpy.sum( numpy.sum(atoms2 * atoms2,axis=0),axis=0)

  #
  # This beautiful step provides the answer.  V and Wt are the orthonormal
  # bases that when multiplied by each other give us the rotation matrix, U.
  # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
  d = numpy.einsum("ijk,jl->ikl", array_atoms1, atoms2)
  V, S, Wt = numpy.linalg.svd( d )

  # We already have our solution, in the results from SVD.
  # we just need to check for reflections and then produce
  # the rotation.  V and Wt are orthonormal, so their det's
  # are +/-1.0 (and thus products are +/- 1.0 ).
  reflect = numpy.linalg.det(V) * numpy.linalg.det(Wt)

  S[:,-1] *= reflect
  V[:,:,-1] *= reflect[:, numpy.newaxis]
  U = numpy.einsum('...ij,...jk->...ki', V, Wt)
  RMSD = E0 - (2.0 * S.sum(axis=1))
  RMSD = numpy.sqrt(abs(RMSD / L))
  return U, COM1-COM2, RMSD

def rmsd(atoms1, atoms2):
  assert len(atoms1) == len(atoms2)
  assert len(atoms1) > 0
  import numpy
  d = atoms1 - atoms2
  d2 = d * d
  sd = d2.sum(axis=1)
  return numpy.sqrt(sd.mean())
