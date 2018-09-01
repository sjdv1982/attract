from __future__ import print_function
import sys, os, tempfile

currdir = os.path.abspath(os.path.split(__file__)[0])

def run_pdb2pqr(pdblines):
  """
  Runs PDB2PQR on the input
  """
  import pdb2pqr
  oldsyspath = sys.path
  pqrhandle, pqrfile = tempfile.mkstemp()
  tmphandle, tmpfile = tempfile.mkstemp()
  tmpf = open(tmpfile, "w")
  for l in pdblines:
    if l.find("XXX") == -1:
      print >> tmpf, l
  tmpf.close()

  try:
    args = [pdb2pqr.__file__, "--ff=charmm", tmpfile, pqrfile]
    #if pdb2pqr.PACKAGE_PATH != "":
    #  sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
    oldstdout = sys.stdout
    sys.stdout = sys.stderr
    oldargv = list(sys.argv)
    sys.argv[:] = args
    pdb2pqr.mainCommand(args)
    sys.argv[:] = oldargv
    sys.stdout = oldstdout
    pqr = os.fdopen(pqrhandle)
    pqrlines = pqr.readlines()
    pqr.close()
  finally:
    os.remove(tmpfile)
    os.remove(pqrfile)

  result = []
  repl = (
    (" H  ", " HN "),
    (" H  ", " HT1"),
    (" H2 ", " HT2"),
    (" H3 ", " HT3"),
  )
  for l in pqrlines:
    result.append(l)
    atom2 = l[12:16]
    for pin, pout in repl:
      if atom2 == pin:
        p = l[:12] + pout + l[16:]
        result.append(p)
  return result

def update_patches(pdb):
  """
  Updates PDB His/Cys patches based on the presence of hydrogens
  """
  for res in pdb:
    if res.resname == "HIS":
      protons = []
      for proton in ("HE1", "HE2", "HD1", "HD2"):
        if proton in res.coords:
          protons.append(proton)

      if protons == ["HE2", "HD1"]: #HIS+ from aareduce
        protons = ["HE1", "HE2", "HD1", "HD2"]
      elif protons == ["HE2"]: #HISE from aareduce
        protons = ["HE1", "HE2", "HD2"]
      elif protons == ["HD1"]: #HISD from aareduce
        protons = ["HE1", "HD1", "HD2"]

      if len(protons) < 3:
        continue
      if "HE2" not in protons or "HE1" not in protons:
        res.topology.patch(top_patches["hisd"])
      if "HD1" not in protons or "HD2" not in protons:
        res.topology.patch(top_patches["hise"])
    if res.resname == "CYS":
      if "HG" not in res.coords:
        res.topology.patch(top_patches["disu"])

def pdbcomplete(pdb, other):
  """
  Completes a PDB by taking missing coordinates from the other PDB
  """
  assert len(pdb) == len(other)
  for r1, r2 in zip(pdb, other):
    for atom in r2.coords:
      if atom not in r1.coords:
        r1.coords[atom] = r2.coords[atom]

def pdbfix(pdb, refe):
  """
  Fixes His and Cys hydrogens in a PDB that are present in the reference
   Hydrogens are filled in by computing the N-H or S-H vector from the reference
   and applying it to the N or S in the reference
  Also copies the OXT position from the reference if it is missing
  """
  assert len(pdb) == len(refe)
  for res, ref in zip(pdb, refe):
    if "oxt" in res.topology.atomorder:
      if "OXT" not in res.coords and "OXT" in ref.coords:
        res.coords["OXT"] = ref.coords["OXT"]
    if res.resname == "HIS":
      for proton in ("HD1", "HE2"):
        if proton in res.coords: continue
        if proton not in ref.coords: continue
        nitro = "N" + proton[1:]
        if nitro not in res.coords or nitro not in ref.coords: #can't fix...
          continue
        nref = ref.coords[nitro]
        href = ref.coords[proton]
        vref = href[0]-nref[0],href[1]-nref[1],href[2]-nref[2]
        nres = res.coords[nitro]
        hres = nres[0] + vref[0], nres[1] + vref[1], nres[2] + vref[2]
        res.coords[proton] = hres
    if res.resname == "CYS":
      proton = "HG"
      sulfur = "SG"
      if proton in res.coords: continue
      if proton not in ref.coords: continue
      if sulfur not in res.coords or sulfur not in ref.coords: #can't fix...
        continue
      sref = ref.coords[sulfur]
      href = ref.coords[proton]
      vref = href[0]-sref[0],href[1]-sref[1],href[2]-sref[2]
      sres = res.coords[sulfur]
      hres = sres[0] + vref[0], sres[1] + vref[1], sres[2] + vref[2]
      res.coords[proton] = hres

class nalib(object):pass

class FixError(Exception):pass

def load_nalib(libname):
  import numpy
  lib = nalib()
  lib.dir = currdir + "/" + libname
  lib.sugar = numpy.load(lib.dir + "/sugar.npy")
  lib.sugaratoms = [l[12:16].strip() for l in open(lib.dir + "/sugar.pdb") if l.startswith("ATOM")]
  lib.phatoms_fit  = ["P", "O1P", "O2P", "O5'", "C5'", "C4'"]  #atoms to fit on
  lib.phatoms_rmsd = ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "C5"]  # to compute best RMSD (avoid base-phosphate clashes)
  lib.phatoms_all  = ["P", "O1P", "O2P", "O5'", "C5'", "C4'"]    # for completion
  lib.base = {}
  lib.baseatoms = {}
  lib.mono = {}
  lib.monoatoms = {}
  lib.ph = {}
  bases = ["A","C","G","U"]
  if libname == "dnalib":
      bases = ["A","C","G","T"]
  for nuc in bases:
    lib.mono[nuc] = numpy.load(lib.dir + "/%s.npy" % nuc)
    lib.monoatoms[nuc] = [l[12:16].strip() for l in open(lib.dir + "/%s.pdb" % nuc) if l.startswith("ATOM")]
    #lib.base[nuc] = numpy.array([(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open(base) if l.startswith("ATOM")])
    lib.base[nuc] = numpy.load(lib.dir + "/base%s.npy" % nuc)
    lib.baseatoms[nuc] = [l[12:16].strip() for l in open(lib.dir + "/base%s.pdb" % nuc) if l.startswith("ATOM")]
    #lib.ph[nuc] = numpy.load(lib.dir + "/%s.npy" % nuc)
    ph = set(lib.phatoms_fit + lib.phatoms_rmsd + lib.phatoms_all)
    lib.ph[nuc] = [ c for a, c in zip(lib.monoatoms[nuc], lib.mono[nuc]) if a in ph]
  return lib

#taken from fit.py
def _apply_matrix(atoms, pivot, rotmat, trans):
  ret = []
  for atom in atoms:
    a = atom-pivot
    atom2 = a.dot(rotmat) + pivot + trans
    ret.append(atom2)
  return ret

def apply_nalib(pdb, lib, heavy, manual):
  """
  Adds missing atoms using a nucleotides lib
  """
  assert heavy == True #TODO: remove this when the rna library has hydrogens
  import numpy
  syspath = list(sys.path)
  sys.path.insert(0, os.environ["ATTRACTTOOLS"])
  import rmsdlib
  from scipy.spatial.distance import cdist
  sys.path[:] = syspath
  for nr, res in enumerate(pdb):
    #print >> sys.stderr, ("residue %i"%nr)
    if res.resname not in ("RA","RC","RG","RU","DA","DC","DG","DT"): continue
    try:
        while 1: #keep fixing as long as we can
            missing = set()
            top = res.topology
            for a in top.atomorder:
                aa = a.upper()
                if heavy and aa[0].startswith("H"): continue
                if aa not in res.coords:
                    missing.add(aa)
            if not missing: break
            #print >> sys.stderr, ("missing", missing)
            nuc = res.resname[1]
            fixmode = None
            for fmode, fatoms in (
                ("ph", lib.phatoms_fit),   #highest priority
                ("sugar", lib.sugaratoms),
                ("base", lib.baseatoms[nuc]),
                ("nucleotide", lib.monoatoms[nuc]),  #lowest priority
                ):
                #we can fix if there are any missing atoms, and there are at least three non-lacking atoms
                fatoms_fit = fatoms
                if fmode == "nucleotide":
                    fatoms_fit = lib.sugaratoms     #in nucleotide mode, fit on the sugar
                if any([(m in fatoms) for m in missing]) and \
                 len([a for a in fatoms_fit if a in res.coords]) >= 3:
                   fixmode = fmode
                   fit_atoms = fatoms_fit
                   break
            if fixmode is None:
                msg = 'residue %s could not be fixed'%res.resid
                raise FixError(msg)
            if fixmode == "base":
                libcoor = lib.base[nuc][numpy.newaxis]
                atoms = lib.baseatoms[nuc]
                #print >> sys.stderr,  lib.baseatoms[nuc], nuc
            elif fixmode == "sugar":
                libcoor = lib.sugar
                atoms = lib.sugaratoms
            elif fixmode == "ph":
                libcoor = lib.ph[nuc]
                atoms = lib.phatoms_all
            elif fixmode == "nucleotide":
                libcoor = lib.mono[nuc]
                atoms = lib.monoatoms[nuc]
            rmsd_atoms = fit_atoms
            #print >> sys.stderr,  fixmode, rmsd_atoms
            #print >> sys.stderr, (fixmode, res.topology)
            if fixmode == "ph":
                rmsd_atoms = lib.phatoms_rmsd
            coor_atoms = [a for a in atoms if a in res.coords and a in fit_atoms]
            coor = numpy.array([res.coords[a] for a in atoms if a in res.coords and a in fit_atoms])
            mask = numpy.array([(a in fit_atoms and a in res.coords) for a in atoms])

            refecoor = numpy.compress(mask, libcoor, axis=1) #or: refecoor = libcoor[:,mask] (slower)
            rotmat, offset, rmsd = rmsdlib.fit(refecoor[0],coor)
            pivot = numpy.sum(coor,axis=0) / float(len(coor))
            fitcoor = numpy.array(_apply_matrix(coor, pivot, rotmat, offset))
            fitcoor = fitcoor.flatten()[numpy.newaxis]
            refecoor = refecoor.reshape((len(refecoor), 3 * len(coor)))

            mask2 = numpy.array([(a in rmsd_atoms and a in res.coords) for a in atoms])
            refecoor2 = numpy.compress(mask2, libcoor, axis=1) #or: refecoor2 = libcoor[:,mask2] (slower)
            coor2 = numpy.array([res.coords[a] for a in atoms if a in res.coords and a in rmsd_atoms])
            rmsdcoor = numpy.array(_apply_matrix(coor2, pivot, rotmat, offset))
            rmsdcoor = rmsdcoor.flatten()[numpy.newaxis]
            refecoor2 = refecoor2.reshape((len(refecoor2), 3 * len(coor)))

            #d = cdist(fitcoor, refecoor, 'sqeuclidean')[0]
            d = cdist(rmsdcoor, refecoor2, 'sqeuclidean')[0]
            libconfnr = numpy.argmin(d)
            libconf = libcoor[libconfnr]
            libconf = _apply_matrix(libconf, pivot+offset, rotmat.T, -offset)
            for anr, a in enumerate(atoms):
                if a in missing or fixmode == "base":
                    #print "FIX", a, fixmode
                    x,y,z = libconf[anr]
                    res.coords[a] = x,y,z
            if fixmode == "nucleotide": break
    except FixError as err:
        if manual:
            e = "\n" + "!"*60 + "\n"
            print(e + "WARNING: " + err.args[0] + e, file=sys.stderr)
        else:
            raise

def pdb_lastresort(pdb):
  """
  Last-resort fixes to prevent errors
  """
  for res in pdb:
    if "oxt" in res.topology.atomorder:
      #Put an OXT at the same place as the O
      if "OXT" not in res.coords and "O" in res.coords:
        res.coords["OXT"] = res.coords["O"]
    if res.resname == "DC":
      #WHATIF refuses to add H3'
      if "H3'" not in res.coords and "C3'" in res.coords:
        res.coords["H3'"] = res.coords["C3'"]
    if res.resname == "DA":
      #WHATIF refuses to add H2, H2', H2'', H1'
      pairs = (("H2","C2"),("H2'","C2'"),("H2''","C2'"),("H1'","C1'"))
      for p1, p2 in pairs:
        if p1 not in res.coords and p2 in res.coords:
          res.coords[p1] = res.coords[p2]
    if res.resname == "RC":
      #WHATIF refuses to add H5, H6
      pairs = (("H5","C5"),("H6","C6"))
      for p1, p2 in pairs:
        if p1 not in res.coords and p2 in res.coords:
          res.coords[p1] = res.coords[p2]
    if "H5T" not in res.coords and "O5'" in res.coords:
      #WHATIF resfuses to add 5TER
      res.coords["H5T"] = res.coords["O5'"]
    if "H3T" not in res.coords and "O3'" in res.coords:
      #WHATIF resfuses to add 3TER
      res.coords["H3T"] = res.coords["O3'"]
