"""
Calculates the (integer) mass of all heavy atoms (C,N,O,S) in a protein file
The protein file can contain one of three things:
- A sequence of amino acids (one-letter code
- A PDB containing ATTRACT atom types (no other forcefield)
- A PDB containing no atom types, but heavy atoms (can also be non-protein, as long as only C,N,O,S is present)
All hydrogen masses are ignored, as well as any isotope distribution
For sequences, the C-terminal oxygen is ignored as well
"""

element_masses = {
  "C": 13.5, #on average, 1.5 hydrogen per C
  "N": 15, #on average, one hydrogen per N
  "O": 16, #most oxygens are backbone (no H)
  "S": 33, #on average, one hydrogen per S
}
skipped = "OXT", "OT2"

atomtype_masses = {
  1: ("CA","C"), #Gly
  2: ("CA","C", "CB"), #Ala
  3: ("CA","C", "CB", "CG"), #Arg
  4: ("CD", "NE","CZ","NH1","NH2"), #Arg
  5: ("CA","C", "CB","CG","OD1","ND2"), #Asn
  6: ("CA","C", "CB","CG","OD1","OD2"), #Asp
  7: ("CA","C", "CB","SG"), #Cys
  8: ("CA","C", "CB", "CG"), #Gln
  9: ("CD", "OE1", "NE2"), #Gln
  10: ("CA","C", "CB", "CG"), #Glu
  11: ("CD", "OE1", "OE2"), #Glu
  12: ("CA","C", "CB", "CG"), #His
  13: ("ND1","CD2", "NE2", "CE1"), #His
  14: ("CA","C", "CB","CG1","CG2","CD1"), #Ile
  15: ("CA","C", "CB","CG","CD1","CD2"), #Leu
  16: ("CA","C", "CB", "CG"), #Lys
  17: ("CD", "CE", "NZ"), #Lys
  18: ("CA","C","CB", "CG"), #Met
  19: ("SD", "CE"), #Met
  20: ("CA","C", "CB", "CG"), #Phe
  21: ("CD1","CD2","CE1","CE2","CZ"), #Phe
  22: ("CA","C", "CB", "CG", "CD"), #Pro
  23: ("CA","C", "CB", "OG"), #Ser
  24: ("CA","C", "CB", "OG1", "CG2"), #Thr
  25: ("CA","C", "CB", "CG"), #Trp
  26: ("CD1","NE1","CD2","CE2","CE3","CH2","CZ3","CZ2"), #Trp
  27: ("CA","C", "CB", "CG"), #Tyr
  28: ("CD1","CD2","CE1","CE2","CZ", "OH"), #Tyr
  29: ("CA","C", "CB", "CG1", "CG2"), #Val
  30: "N",
  31: "O",
}

aa_masses = {
  "A": ("CA","C","O","N","CB"),
  "C": ("CA","C","O","N","CB", "SG"),
  "D": ("CA","C","O","N","CB", "CG","OD1","OD2"),
  "E": ("CA","C","O","N","CB", "CG", "CD", "OE1", "OE2"),
  "F": ("CA","C","O","N","CB", "CG", "CD1","CD2","CE1","CE2","CZ"),
  "G": ("CA","C","O","N"),
  "H": ("CA","C","O","N","CB", "CG", "ND1","CD2", "NE2", "CE1"),
  "I": ("CA","C","O","N","CB", "CG1","CG2","CD1"),
  "K": ("CA","C","O","N","CB", "CG", "CD", "CE", "NZ"),
  "L": ("CA","C","O","N","CB", "CG","CD1","CD2"),
  "M": ("CA","C","O","N","CB", "CG", "SD", "CE"),
  "N": ("CA","C","O","N","CB", "CG","OD1","ND2"),
  "P": ("CA","C","O","N","CB", "CG", "CD"),
  "Q": ("CA","C","O","N","CB", "CG", "CD", "OE1", "NE2"),
  "R": ("CA","C","O","N","CB", "CG", "CD", "NE","CZ","NH1","NH2"),
  "S": ("CA","C","O","N","CB", "OG"),
  "T": ("CA","C","O","N","CB", "OG1", "CG2"),
  "V": ("CA","C","O","N","CB", "CG1", "CG2"),
  "W": ("CA","C","O","N","CB", "CG", "CD1","NE1","CD2","CE2","CE3","CH2","CZ3","CZ2"),
  "Y": ("CA","C","O","N","CB", "CG", "CD1","CD2","CE1","CE2","CZ", "OH"),
}
num = "0","1","2","3","4","5","6","7","8","9"

import sys
proteinfile = sys.argv[1]
lines = open(proteinfile).readlines()

#1. Are we an amino acid sequence?
ok = True
sequence = ""
for l in lines:
  if not ok: break
  if l.lstrip().startswith(">"): continue
  for ll in l:
    if len(ll.strip()) == 0: continue
    elif ll in num: continue
    elif ll in aa_masses: sequence = sequence + ll
    else:
      ok = False
      break

if ok: #file contains amino acid sequence
  mass = 0
  for aa in sequence:
    for atom in aa_masses[aa]:
      mass += element_masses[atom[0]]
  print(int(mass))
  sys.exit()
      
#2. Are we a PDB with ATTRACT atom types?
ok = True
mass = 0
for l in lines:  
  if not l.startswith("ATOM"): continue
  if len(l) < 59:
    ok = False
    break
  atomtype = l[57:59]
  try:
    atomtype = int(atomtype)
  except ValueError:
    ok = False
    break  
  if atomtype in (32,99): continue  
  if atomtype not in atomtype_masses:
    if l[56] == " " and l[59] == " ":
      raise Exception("PDB file has not been reduced using ATTRACT forcefield for proteins")
    else:
      ok = False
      break
  
  atoms = atomtype_masses[atomtype]
  if isinstance(atoms, str): atoms = [atoms]
  for a in atoms:
    mass += element_masses[a[0]]
if ok: #file contains reduced PDB
  print(int(mass))
  sys.exit()
  

for l in lines:  
  if not l.startswith("ATOM"): continue
  if len(l) < 55: raise Exception("Malformed PDB")
  if l[12] == "H": continue
  code = l[13:16]
  if code in skipped: continue
  element = code[0]
  if element == " ": raise Exception("Malformed PDB")
  if element not in element_masses: raise Exception("Unknown element '%s'" % element)
  mass += element_masses[element]
print(int(mass))