def determine_passive_residues(activereslist, surfacelist, pdb, radius):
  passivereslist0 = set()  
  distances = pdb.convert(InterResidueDistanceArray)    
  for d in distances:
    if d.distance > radius: continue    
    if d.residue1 in activereslist:
      if d.residue2 in surfacelist: passivereslist0.add(d.residue2)
    elif d.residue2 in activereslist:
      if d.residue1 in surfacelist: passivereslist0.add(d.residue1)
  passivereslist = []
  for r in passivereslist0:
    if r not in activereslist: passivereslist.append(r)  
  passivereslist.sort()
  return passivereslist
