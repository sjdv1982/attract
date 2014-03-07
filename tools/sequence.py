import sys
pdb = sys.argv[1]

code = {
  "ALA": "A",
  "CYS": "C",
  "ASP": "D",
  "GLU": "E",
  "PHE": "F",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LYS": "K",
  "LEU": "L",
  "MET": "M",
  "ASN": "N",
  "PRO": "P",
  "GLN": "Q",
  "ARG": "R",
  "SER": "S",
  "THR": "T",
  "VAL": "V",
  "TRP": "W",
  "TYR": "Y",
}

lines = open(pdb).readlines()
id = None
seq = ""
for l in lines:
  if not l.startswith("ATOM"): continue
  cid = l[17:26]
  if id != cid:
    id = cid
    aa = l[17:20]
    seq += code[aa]

print(seq)    