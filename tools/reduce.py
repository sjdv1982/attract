import sys, os
reducedat = os.path.split(os.path.abspath(__file__))[0] + os.sep + "reduce.dat"

def read_forcefield(forcefieldfile):
  ff = {}
  aa = None
  for l in open(forcefieldfile):
    pound = l.find("#")
    if pound > -1: l = l[:pound]
    l = l.strip()
    if not len(l): continue
    ll = l.split()
    if len(ll) == 1:
      aa = ll[0]
      assert len(aa) <= 3, l
      ff[aa] = []
    else:
      assert aa is not None
      try:
        atomtype = int(ll[0])
      except ValueError:
        raise ValueError(l)
      atoms = ll[2:]
      charge = 0.0
      try:
        charge = float(atoms[-1])
        atoms = atoms[:-1]
      except ValueError:
        pass
      ff[aa].append( (int(ll[0]), ll[1], set(atoms), charge) )
  return ff  
      
  
ff = read_forcefield(reducedat)
prot = sys.argv[1]

res = None
resname = None
rescounter = 0
atomcounter = 0
rescoor = {}

def print_res():
  global rescounter, atomcounter, rescoor
  if not len(rescoor): return  
  rescounter += 1
  for l in ff[resname]:
    if (l[0], l[1]) not in rescoor: continue
    c = rescoor[(l[0], l[1])]
    x, y, z = c[1]/c[0], c[2]/c[0], c[3]/c[0]
    """
    KLUDGES to reproduce old reduce behavior:    
    """
    if resname in ("GLU", "GLN") and l[1] in ("CN1", "CO1"):
      x, y, z = c[1]*0.333, c[2]*0.333, c[3]*0.333
    """
    /KLUDGES
    """
    atomcounter += 1
    line = (atomcounter, l[1], resname, rescounter, x, y, z, l[0], l[3], 1.0)
    print "ATOM%7d  %-3s %-3s  %4d    %8.3f%8.3f%8.3f%5d%8.3f 0%5.2f" % line
  rescoor = {}
  
for l in open(prot):
  if not l.startswith("ATOM"): continue
  cres = l[21:26]
  if cres != res:
    print_res()
    res = cres
    resname = l[17:20].strip()
    assert resname in ff, l
    ffres = ff[resname] 
  try:  
    atom = l[12:16].strip()
    x = float(l[30:38])
    y = float(l[38:46])
    z = float(l[46:54])
  except ValueError:
    continue
  for bead in ffres:
    for at in bead[2]:
      if atom != at: continue
      beadname = bead[0], bead[1]
      if beadname not in rescoor: 
        rescoor[beadname] = [0, 0.0, 0.0, 0.0]
      c = rescoor[beadname]
      c[0] += 1
      c[1] += x
      c[2] += y
      c[3] += z
      break
print_res()    