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
  "C": 13.5, #on average, 1.5 hydrogen per C
  "N": 15, #on average, one hydrogen per N
  "O": 16, #most oxygens are backbone (no H)
  "S": 33, #on average, one hydrogen per S
}


saxs_factors = { #vacuum values from imp/modules/saxs/src/FormFactorTable.cpp 
  "C": 4.99 + 1.5 * 0.999953,  #on average, 1.5 hydrogen per C
  "N": 5.9992 + 1.0 * 0.999953, #on average, one hydrogen per N
  "O": 6.9946, #most oxygens are backbone (no H)
  "S": 14.9993 + 1.0 * 0.999953, #on average, one hydrogen per S
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
      w = element_masses[e]
    weights.append(w)
  return np.array(weights,dtype="float32")

def get_saxs_factors(pdb):  
  factors = []
  pdb = rmsdlib.read_pdb(pdb)
  for a in pdb.atoms():
    atomtype = int(a.line[57:59])
    w = 0.0
    if a.name in skipped: continue
    e = a.name[0]
    if e == "H": continue
    f = saxs_factors[e]
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
        if e == "H": continue
        w = element_masses[e]
      weights.append(w)
    arr.append(np.array(weights,dtype="float32"))
  return arr  

def get_coordinates(pdb, ignore_weightless):
  pdb = rmsdlib.read_pdb(pdb)
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
  