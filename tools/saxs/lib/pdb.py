import sys, os
sys.path.append(os.environ["ATTRACTTOOLS"])
import rmsdlib
import numpy as np

aweights = (
 0, #0
 27.0, #1
 40.5, #2
 54.0, #3
 72.0, #4
 85.0, #5
 86.0, #6
 73.5, #7
 54.0, #8
 44.5, #9
 54.0, #10
 45.5, #11
 54.0, #12
 57.0, #13
 81.0, #14
 81.0, #15
 54.0, #16
 42.0, #17
 54.0, #18
 46.5, #19
 54.0, #20
 67.5, #21
 67.5, #22
 56.5, #23
 70.0, #24
 54.0, #25
 109.5, #26
 54.0, #27
 83.5, #28
 67.5, #29
 15, #30
 16, #31
)

element_masses = {
  "H": 1.0,
  "C": 12, #on average, 1.5 hydrogen per C
  "N": 14, #on average, one hydrogen per N
  "O": 16, #most oxygens are backbone (no H)
  "S": 32, #on average, one hydrogen per S
}

def get_hydrogens(resname, name):
  #Adapted from IMP saxs module modules/saxs/src/FormFactorTable.cpp
  if name == 'N':
    if resname == 'PRO':
      return 0.0
    else:
      return 1.0
  elif name == 'CA':
    if resname == 'GLY':
      return 2.0  
    else:
      return 1.0
    
      
  elif name == 'C':
    return 0.0
  elif name == 'CH2':
    return 2.0
  
  elif name == 'CH3':
    return 3.0
  
  elif name == 'CH':
    return 1.0
  
  elif name == 'OH':
    return 1.0
  
  elif name == 'NH':
    return 1.0
  
  elif name == 'NH2':
    return 2.0
  
  elif name == 'CB':
    if resname in ['ILE','THR','VAL']:
      return 1.0
    elif resname == 'ALA':
      return 3.0
    else:
      return 2.0
    
  elif name == 'CG':
    if resname in ['ASN','HIS','ASP','PHE','TRP','TYR']:
      return 0.0
    elif resname == 'LEU':
      return 1.0
    else:
      return 2.0
    
  elif name == 'CG1':
    if resname == 'ILE':
      return 2.0
    elif resname == 'VAL':
      return 3.0
    
  elif name == 'CG2':
    return 3.0
  
  elif name == 'CD':
    if resname in ['GLU','GLN']:
      return 0.0
    else:
      return 2.0
    
  elif name == 'CD1':
    if resname in ['ILE','LEU']:
      return 3.0
    elif resname in ['PHE','TRP','TYR']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'CD2':
    if resname == 'LEU':
      return 3.0
    elif resname in ['PHE','HIS','TYR']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'CE':
    if resname == 'LYS':
      return 2.0
    elif resname == 'MET':
      return 3.0
    else:
      return 0.0
    
  elif name == 'CE1':
    if resname in ['PHE','HIS','TYR']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'CE2':
    if resname in ['PHE','TYR']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'CZ':
    if resname == 'TRP':
      return 1.0
    else:
      return 0.0
    
  elif name in ['CZ1','CZ2','CZ3','CE3']:
    if resname == 'TRP':
      return 1.0
    else:
      return 0.0
    
  #DNA, TODO  check naming, check
  elif name == "C5'":
    return 2.0
  
  elif name in ["C1'","C2'","C3'","C4'"]:
    return 1.0
  
  elif name == 'C2':
    if resname in ['DA','RA']:
     return 1.0
    else:
      return 0.0
    
  elif name == 'C4':
    return 0.0
  
  elif name == 'C5':
    if resname in ['DC','RC','DU','RU']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'C6':
    if resname in ['DC','RC','DU','RU','DT','RT']:
      return 1.0
    else:
      return 0.0
    
  elif name == 'C7':
    return 3.0
  elif name == 'C8':
    return 1.0
  
  elif name == 'ND1':
    if resname == 'HIS':
      return 1.0
    else:
      return 0.0
    
  elif name == 'ND2':
    if resname == 'ASN':
      return 2.0
    else:
      return 0.0
    
  elif name == 'NH1' or name == 'NH2':
    if resname == 'ARG':
      return 2.0
    else:
      return 0.0
  
  elif name == 'NE':
    if resname == 'ARG':
      return 1.0
    else:
      return 0.0
    
  elif name == 'NE1':
    if resname == 'TRP':
      return 1.0
    else:
      return 0.0
    
  elif name == 'NE2':
    if resname == 'GLN':
      return 2.0
    else:
      return 0.0
    
  elif name == 'NZ':
    if resname == 'LYS':
      return 3.0
    else:
      return 0.0
    
  elif name == 'N1':
    if resname in ['DG','RG']:
      return 1.0
    else:
      return 0.0
  
  elif name in ['N2','N4','N6']:
    return 2.0
  
  elif name == 'N3':
    if resname in ['DU','RU']:
      return 1.0
    else:
      return 0.0
    
  elif name in ['N7','N9']:
    return 0.0
  
  elif name in ['O','OE1','OE2','OD1','OD2','O1A','O2A','OXT','OT1','OT2']:
    return 0.0
  
  elif name == 'OG':
    if resname == 'SER':
      return 1.0
    else:
      return 0.0
    
  elif name == 'OG1':
    if resname == 'THR':
      return 1.0
    else:
      return 0.0
    
  elif name == 'OH':
    if resname == 'TYR':
      return 1.0
    else:
      return 0.0
    
  elif name in ['O1P',"O3'","O2P","O4'","O5'","O2","O4","O6"]:
    return 0.0
  
  elif name == "O2'":
    return 1.0
  
  elif name == 'HOH':
    return 2.0
  
  elif name == 'SD':
    return 0.0
  
  elif name == 'SG':
    if resname == 'CYS':
      return 1.0
    else:
      return 0.0
    
  else:
    print "ERROR Atom name", name, "in residue",resname,"not recognized!"
    sys.exit(1)
      
  
saxs_factors = { #vacuum-dummy values from imp/modules/saxs/src/FormFactorTable.cpp 
  "H": [0.999953,1.7201],
  "C": [5.9992,5.49096],  #on average, 1.5 hydrogen per C
  "N": [6.9946,0.83166], #on average, one hydrogen per N
  "O": [7.9994,3.04942], #most oxygens are backbone (no H)
  "S": [15.9998,6.63324], #on average, one hydrogen per S
}

skipped = "OXT", "OT2"

def get_weights(pdb, ignore_weightless):  
  weights = []
  pdb = rmsdlib.read_pdb(pdb)
  for a in pdb.atoms():
    atomtype = int(a.line[57:59])
    w = 0.0
    if atomtype in (32,99) and ignore_weightless:
      continue
    elif atomtype < 32:
      w = aweights[atomtype]
    elif atomtype in (32,99):
      pass
    else:
      if a.name in skipped: continue
      e = a.name[0]
      if e == "H": continue
      hnum = get_hydrogens(a.resname,a.name)
      w = element_masses[e]+hnum*element_masses["H"]
    weights.append(w)
  return np.array(weights,dtype="float32")

def get_vacuum_factors(pdb):  
  factors = []
  pdb = rmsdlib.read_pdb(pdb)
  for a in pdb.atoms():
    atomtype = int(a.line[57:59])
    w = 0.0
    if a.name in skipped: continue
    e = a.name[0]
    if e == "H": continue
    hnum = get_hydrogens(a.resname,a.name)
    f = saxs_factors[e][0]+hnum*saxs_factors["H"][0] #vacuum form factors
    factors.append(f)
  return np.array(factors,dtype="float32")

def get_dummy_factors(pdb):  
  factors = []
  pdb = rmsdlib.read_pdb(pdb)
  for a in pdb.atoms():
    atomtype = int(a.line[57:59])
    w = 0.0
    if a.name in skipped: continue
    e = a.name[0]
    if e == "H": continue
    hnum = get_hydrogens(a.resname,a.name)
    f = saxs_factors[e][1]+hnum*saxs_factors["H"][1] #dummy solvent form factors
    factors.append(f)
  return np.array(factors,dtype="float32")

def get_weights_multi_attract(pdb, ignore_weightless):  
  pdbs = rmsdlib.read_multi_attract_pdb(pdb)
  arr = []
  for pdb in pdbs:
    weights = []
    for a in pdb.atoms():
      atomtype = int(a.line[57:59])
      w = 0.0
      if atomtype in (32,99) and ignore_weightless:
        continue
      elif atomtype < 32:
        w = aweights[atomtype]
      elif atomtype in (32,99):
        pass
      else:
        if a.name in skipped: continue
        e = a.name[0]
        hnum = hydrogens[a.name]
        if a.resname == 'GLY' and a.name == 'CA':
          hnum = 2.0
        if a.resname == 'PRO' and a.name == 'N':
          hnum = 0.0
        if e == "H": continue
        w = element_masses[e]+hnum*element_masses["H"]
      weights.append(w)
    arr.append(np.array(weights,dtype="float32"))
  return arr  

def get_coordinates(pdb, ignore_weightless):
  pdb = rmsdlib.read_pdb(pdb)
  if ignore_weightless:
    coors = []  
    for a in pdb.atoms():
      if a.name in skipped: continue
      e = a.name[0]
      if e == "H": continue
      coors.append((a.x,a.y,a.z))    
  else:
    coors = pdb.coordinates()
  return np.array(coors, dtype="float32")
    
def get_coordinates_multi_attract(pdb, ignore_weightless):
  pdbs = rmsdlib.read_multi_attract_pdb(pdb)
  arr = []
  for pdb in pdbs:
    if ignore_weightless:
      coors = []  
      for a in pdb.atoms():
        atomtype = int(a.line[57:59])
        if atomtype in (32,99):
          continue
        elif atomtype < 32:
          pass
        else:
          if a.name in skipped: continue
          e = a.name[0]
          if e == "H": continue
        coors.append((a.x,a.y,a.z))    
    else:
      coors = pdb.coordinates()
    arr.append(np.array(coors, dtype="float32"))
  return arr  
  
