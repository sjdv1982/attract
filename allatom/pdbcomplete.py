#!/usr/bin/env python2

from __future__ import print_function
import sys, os, tempfile, numpy as np

currdir = os.path.abspath(os.path.split(__file__)[0])

na_resnames = {"RA","RC","RG","RU","DA","DC","DG","DT","A","T","C","G","U"}
# pdb2pqr returns another atom name for that RNA atom than the opls ff.
# check_pdb will chock on it.
map_atnames = {"HO2'": "H2''"}

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
      print(l, file=tmpf)
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

def update_patches(pdb, top_patches):
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
        res.topology.patch(top_patches["HISD"])
      if "HD1" not in protons or "HD2" not in protons:
        res.topology.patch(top_patches["HISE"])
    if res.resname == "CYS":
      if "HG" not in res.coords:
        res.topology.patch(top_patches["DISU"])

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
    lib = nalib()
    lib.dir = currdir + "/" + libname
    lib.sugar = {}
    lib.base = {}
    lib.nucl = {}
    lib.ph = {}
    lib.pho5 = {}
    ph = {
        # !!! keep the order of atoms in the library!
        "atoms":     ["P", "O1P", "O2P", "O5'"], # for completion
        #"fit_atoms": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "C3'"], #atoms to fit on
        "fit_atoms": ["P", "O1P", "O2P", "O5'"], #atoms to fit on
        "rmsd_atoms": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "N3", "C5", "C3'"], # to compute best RMSD (avoid base-phosphate clashes),
        "max_missing": 4,
        }
    pho5 = {
    "coor": np.load(lib.dir + "/5PHO.npy"),
    "atoms": [l[12:16].strip() for l in open(lib.dir + "/5PHO.pdb") if l.startswith("ATOM")],
    "max_missing": 5,
    }
    pho5["fit_atoms"] = pho5["atoms"]
    pho5["rmsd_atoms"] = pho5["atoms"]
    sugar = {
        "coor": np.load(lib.dir + "/sugar.npy"),
        "atoms": [l[12:16].strip() for l in open(lib.dir + "/sugar.pdb") if l.startswith("ATOM")],
        "max_missing": 2,
        }
    sugar["fit_atoms"] = sugar["atoms"]
    sugar["rmsd_atoms"] = sugar["atoms"]
    bases = ["A","C","G","U"]
    if libname == "dnalib":
        bases = ["A","C","G","T"]
    for nuc in bases:
        lib.sugar[nuc] = sugar
        lib.nucl[nuc] = {
            "coor": np.load(lib.dir + "/%s.npy" % nuc),
            "missing": 23,
        }
        nuclatoms = [l[12:16].strip() for l in open(lib.dir + "/%s.pdb" % nuc) if l.startswith("ATOM")]
        lib.nucl[nuc]["atoms"] = nuclatoms
        lib.nucl[nuc]["fit_atoms"] = ["C4'", "O4'", "C1'", "C2'", "C3'"]
        lib.nucl[nuc]["rmsd_atoms"] = nuclatoms

        lib.base[nuc] = {}
        #lib.base[nuc]["coor"] = np.array([(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open(base) if l.startswith("ATOM")])
        lib.base[nuc]["coor"] = np.load(lib.dir + "/base%s.npy" % nuc)[None]
        baseatoms = [l[12:16].strip() for l in open(lib.dir + "/base%s.pdb" % nuc) if l.startswith("ATOM")]
        lib.base[nuc]["atoms"] = baseatoms
        lib.base[nuc]["fit_atoms"] = baseatoms
        lib.base[nuc]["rmsd_atoms"] = baseatoms
        lib.base[nuc]["max_missing"] = len(baseatoms) - 3

        lib.ph[nuc] = ph.copy()
        lib.ph[nuc]["coor"] = np.load(lib.dir + "/%s.npy" % nuc)
        phlist = set(lib.ph[nuc]["atoms"] + lib.ph[nuc]["rmsd_atoms"] + lib.ph[nuc]["fit_atoms"])
        ph_indices = [ anr for anr,a in enumerate(lib.nucl[nuc]["atoms"]) if a in phlist]
        #lib.ph[nuc]["coor"] = lib.nucl[nuc]["coor"][:, ph_indices]

        lib.pho5[nuc] = pho5
    return lib

def check_missing(pdb):
    for nr, res in enumerate(pdb):
        top = res.topology
        for a in top.atomorder:
            aa = a.upper()
            if aa[0].startswith("H"): continue
            if aa not in res.coords:
                return 1
    return 0

def process_pdb(pdb):
    other_coor, atom_to_residue = [], []
    for nr, res in enumerate(pdb):
        top = res.topology
        for a in top.atomorder:
            aa = a.upper()
            if aa[0].startswith("H"): continue
            if aa not in res.coords:
                continue
            xyz = res.coords[aa]
            other_coor.append(xyz)
            atom_to_residue.append(nr)
    return np.array(other_coor), np.array(atom_to_residue)

#taken from fit.py
def _apply_matrix(atoms, pivot, rotmat, trans):
  ret = []
  for atom in atoms:
    a = atom-pivot
    atom2 = a.dot(rotmat) + pivot + trans
    ret.append(atom2)
  return ret

def _apply_matrix_multi(atoms_array, pivots, rotmats, offsets):
    cen = (atoms_array - pivots[:, None, :])
    newcoor = np.einsum("ijk,ikl->ijl", cen, rotmats) #diagonally broadcasted form of cen.dot(rotmats)
    newcoor += (pivots + offsets)[:, None, :]
    return newcoor

def check_clashes(nl, libconf, lib_complete_indices, new_at, tree, atom_to_residue):
    #for nl, libconf in enumerate(libcoor_fitted_sorted):
    lib_complete = libconf[lib_complete_indices]
    # get clashes with any original residue
    neighbors = tree.query_ball_point(lib_complete, r=2)
    neighbors = np.concatenate(neighbors).astype(int)
    # ! the completed residue clashes w. himself
    clash_res = np.unique(atom_to_residue[neighbors])
    # get clashes with other added atom
    new_clash = 0
    if len(new_at):
        new_atoms = np.array(new_at)
        dd = new_atoms[None,:,:] - lib_complete[:,None,:]
        dist_new = np.sum(dd*dd, axis=2)
        new_clash = np.sum(dist_new < 4)
        #m = np.min(dist_new)
        #print('nb new clash: %i, min %f'%(new_clash, m))
    if len(clash_res) < 2 and not new_clash:
        return 0
    return 1

def rank_library(res, libcoor, fit_atoms, rmsd_atoms, atoms):
    '''
    rank nucleotides in library by best-fitting to our target residue
    '''
    import rmsdlib
    from scipy.spatial.distance import cdist
    #select atoms in the order of the library atoms
    coor_atoms = [a for a in atoms if a in res.coords and a in fit_atoms]
    coor = np.array([res.coords[a] for a in atoms if a in res.coords and a in fit_atoms])
    fit_mask = np.array([(a in fit_atoms and a in res.coords) for a in atoms])
    #fit the mononucl lib on the nucl to repair:
    ## get the rotation-translation to apply to the rank_library
    libcoor_fit = libcoor[:,fit_mask]
    rotmats, offsets, rmsds = rmsdlib.multifit(libcoor_fit,coor)
    rotmats = rotmats.swapaxes(1,2)
    offsets = -offsets
    x = rmsds.argsort()
    pivots = libcoor_fit.sum(axis=1)/libcoor_fit.shape[1]
    ## fit on fit_atoms and evaluate the rmsd on the rmsd_atoms
    common_atoms = set([a for a in rmsd_atoms if a in res.coords])
    rmsd_mask = np.array([(a in common_atoms) for a in rmsd_atoms])
    libcoor_rmsd_unfitted = libcoor[:,rmsd_mask]
    libcoor_rmsd_fitted = _apply_matrix_multi(libcoor_rmsd_unfitted, pivots, rotmats, offsets)
    coor_rmsd = np.array([res.coords[a] for a in rmsd_atoms if a in common_atoms])
    dist = libcoor_rmsd_fitted - coor_rmsd
    d = np.einsum("ijk,ijk->i", dist, dist)
    ### d = rmsds ### use the fit rmsds as an approx. for debugging
    # sort library indices by RMSD
    lib_indices = np.argsort(d)
    libcoor_fitted = _apply_matrix_multi(libcoor, pivots, rotmats, offsets)
    libcoor_fitted_sorted = libcoor_fitted[lib_indices]
    return libcoor_fitted_sorted

def apply_nalib(pdb, lib, manual, heavy=True):
    """
    Adds missing atoms using a nucleotides lib
    """
    print("apply_nalib", file=sys.stderr)
    #assert heavy == True #TODO: remove this when the rna library has hydrogens
    syspath = list(sys.path)
    if "ATTRACTTOOLS" in os.environ:
        sys.path.insert(0, os.environ["ATTRACTTOOLS"])
    import rmsdlib
    from scipy.spatial import cKDTree
    sys.path[:] = syspath
    any_missing = False
    #check if any missing atoms in the whole PDB
    #If any atom missing, then process the PDB to get all other atoms
    any_missing = check_missing(pdb)
    if any_missing:
        #coordinates toward which one should test clashes
        other_coor, atom_to_residue = process_pdb(pdb)
        tree  = cKDTree(other_coor)
    new_at_all = []
    for nr, res in enumerate(pdb):
        if res.resname not in na_resnames: continue
        new_at = []
        try:
            while 1: #keep fixing as long as we can
                missing = set()
                top = res.topology
                for a in top.atomorder:
                    aa = a.upper()
                    #if aa[0].startswith("H"): continue ###
                    if aa.startswith("H") or aa == '5pho': continue ###
                    if aa not in res.coords:
                        missing.add(aa)
                if not missing: break
                nuc = res.resname[1]
                fixmode = None
                for fixmode in ("ph", "sugar", "base", "nucl"): #from high to low priority
                    #we can fix if there are any missing atoms, and there are at least three non-lacking atoms
                    #if fixmode == "sugar": continue
                    sublib = getattr(lib, fixmode) # lib.ph or lib.sugar or ...
                    atoms = sublib[nuc]["atoms"]
                    fit_atoms = sublib[nuc]["fit_atoms"]
                    rmsd_atoms = sublib[nuc]["rmsd_atoms"]
                    libcoor = sublib[nuc]["coor"]
                    if any([(m in atoms) for m in missing]) and \
                     len([a for a in fit_atoms if a in res.coords]) >= 3:
                        break
                else:
                    msg = 'residue %s could not be fixed'%res.resid
                    print("missing atoms: ", file=sys.stderr)
                    print(missing, file=sys.stderr)
                    break
                    #raise FixError(msg)
                print("fixing %s"%fixmode, file=sys.stderr)
                libcoor_fitted_sorted = rank_library(res, libcoor, fit_atoms, rmsd_atoms, atoms)
                lib_complete_indices = []
                for anr, a in enumerate(atoms):
                    if a in missing or fixmode == "base":
                        lib_complete_indices.append(anr)
                #TODO: change clashing threshold when not --heavy
                #optimize: if clashes, take next nucleotide in mononucl_library
                print('optimize resid %s %s'%(res.resname, res.resid), file=sys.stderr)
                for nl, libconf in enumerate(libcoor_fitted_sorted):
                    clashing = check_clashes(nl, libconf, lib_complete_indices, new_at_all, tree, atom_to_residue)
                    if not clashing:
                        break
                if nl == len(libconf):
                    raise FixError("all lib conformer clash on resid %s %s"%(res.resname, res.resid))
                #
                for anr, a in enumerate(atoms):
                    if anr in lib_complete_indices: # add "or 1" to replace all atoms, for debugging
                        x,y,z = libconf[anr]
                        res.coords[a] = x,y,z
                        new_at.append([x,y,z])
                if fixmode == "nucleotide" : break
        except FixError as err:
            if manual:
                e = "\n" + "!"*60 + "\n"
                print(e + "WARNING: " + err.args[0] + e, file=sys.stderr)
            else:
                raise
        for a in new_at:
            new_at_all.append(a)

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
