import sys, os, tempfile
import parse_cns_top

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

def run_whatif(pdblines):
  """
  Runs the WHATIF server on the input
  """  
  from whatif import run  
  whatifresult0 = run("corall", "\n".join(list(pdblines)))
  whatifresult = run("htopo", whatifresult0)
  whatiflines = whatifresult.split("\n")  
  result = []
  repl = (
    (" H  ", " HN "),
    (" H  ", " HT1"),
    (" H2 ", " HT2"),
    (" H3 ", " HT3"),
    ("1HD2", "HD21"),
    ("2HD2", "HD22"),    
    ("1HE2", "HE21"),
    ("2HE2", "HE22"),
    ("1HH1", "HH11"),
    ("2HH1", "HH12"),        
    ("1HH2", "HH21"),
    ("2HH2", "HH22"),
    ("'HO2", "HO2'"),
    (" HO2", "HO2'"),
    ("H2''", " H2'"),    
  )  
  for l in whatiflines:
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
        res.topology.patch(parse_cns_top.presidues["hisd"])
      if "HD1" not in protons or "HD2" not in protons: 
        res.topology.patch(parse_cns_top.presidues["hise"])
    if res.resname == "CYS":
      if "HG" not in res.coords:
        res.topology.patch(parse_cns_top.presidues["disu"])

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
      #WHATIF refuses to add H2, H2', H2''
      pairs = (("H2","C2"),("H2'","C2'"),("H2''","C2'"))
      for p1, p2 in pairs:
        if p1 not in res.coords and p2 in res.coords:
          res.coords[p1] = res.coords[p2]
    if res.resname == "RC":
      #WHATIF refuses to add H5, H6
      pairs = (("H5","C5"),("H6","C6"))
      for p1, p2 in pairs:
        if p1 not in res.coords and p2 in res.coords:
          res.coords[p1] = res.coords[p2]
        
      
