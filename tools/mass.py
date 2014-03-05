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
  "C": 12,
  "N": 14,
  "O": 16,
  "S": 32,
}
skipped = "OXT", "OT2"

atomtype_masses = {
  1: 0,
  2: 0,
  3: 0,
  4: 0,
  5: 0,
  6: 0,
  7: 0,
  8: 0,
  9: 0,
  10: 0,
  11: 0,
  12: 0,
  13: 0,
  14: 0,
  15: 0,
  16: 0,
  17: 0,
  18: 0,
  19: 0,
  20: 0,
  21: 0,
  22: 0,
  23: 0,
  24: 0,
  25: 0,
  26: 0,
  27: 0,
  28: 0,
  29: 0,
  30: "N",
  31: "O",
}

aa_masses = {
  "A": 0,
  "C": 0,
  "D": 0,
  "E": 0,
  "F": 0,
  "G": 0,
  "H": 0,
  "I": 0,
  "K": 0,
  "L": 0,
  "M": 0,
  "N": 0,
  "P": 0,
  "Q": 0,
  "R": 0,
  "S": 0,
  "T": 0,
  "V": 0,
  "W": 0,
  "Y": 0,
}
num = "0","1","2","3","4","5","6","7","8","9"

import sys
proteinfile = sys.argv[1]
lines = open(proteinfile)

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
  for aa in seq:
    mass += aa_masses[aa]
  print(mass)  
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
  for a in atoms:
    mass += element_masses[a]
if ok: #file contains reduced PDB
  print(mass)
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
print(mass)