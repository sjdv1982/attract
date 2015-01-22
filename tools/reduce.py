"""
Reduce converts a PDB into a modified PDB format, where the atom type for each atom is indicated
This reduced PDB is what is understood by ATTRACT.
These atoms are actually pseudo-atoms (beads), representing the average position of one or more real atoms.
The pseudo-atoms definitions are in reduce.dat, see that file for more details.

DNA bases are DA, DC, DG, DT; and RNA bases are RA, RC, RG, RU. One-letter and three-letter nucleic acid\
codes can be  automatically interpreted as DNA or RNA with the --dna and --rna options.

Author: Sjoerd de Vries, Technische Universitaet Muenchen
"""

import sys, os

has_argparse = False
try:
  import argparse  
  has_argparse = True  
except ImportError:
  import optparse  #Python 2.6


#Mapping of nucleic-acid codes to DNA/RNA
mapnuc = {
  "A": ["DA", "RA"],
  "ADE": ["DA", "RA"],
  "C": ["DC", "RC"],
  "CYT": ["DC", "RC"],
  "G": ["DG", "RG"],
  "GUA": ["DG", "RG"],
  "T": ["DT", None],
  "THY": ["DT", None],
  "U": [None, "RU"],
  "URA": [None, "RU"],
  "URI": [None, "RU"],  
} 

if has_argparse:
  parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("pdb",help="PDB file to reduce")
  parser.add_argument("output",help="reduced output PDB file", nargs="?")
else:
  parser = optparse.OptionParser()
  parser.add_argument = parser.add_option
parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
parser.add_argument("--compat",help="Maximize compatibility with reduce.f", action="store_true")

if has_argparse:
  args = parser.parse_args()
else:
  args, positional_args = parser.parse_args()
  args.pdb = None
  args.output = None
  if positional_args:
    if len(positional_args) > 0: args.pdb = positional_args[0]
    if len(positional_args) > 1: args.output = positional_args[1]

if args.dna and args.rna:
  raise ValueError("Options --dna and --rna are mutually exclusive")

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
pdb = args.pdb
assert os.path.exists(pdb)
if args.output is None:
  args.output = os.path.splitext(pdb)[0] + "r.pdb"
outp = open(args.output, "w")  


res = None
resname = None
rescounter = 0
atomcounter = 0
rescoor = {}

def print_res():
  global rescounter, atomcounter, rescoor
  if not len(rescoor): return  
  rescounter += 1
  #print res[1:].strip(), rescounter
  for l in ff[resname]:
    if (l[0], l[1]) not in rescoor: continue
    c = rescoor[(l[0], l[1])]
    x, y, z = c[1]/c[0], c[2]/c[0], c[3]/c[0]
    if args.compat and resname in ("GLN", "GLU") and l[1] in ("CN1", "CO1"):
      x, y, z = c[1] * 0.333, c[2] * 0.333, c[3] * 0.333
    atomcounter += 1
    line = (atomcounter, l[1], resname, rescounter, x, y, z, l[0], l[3], 1.0)
    print >> outp, "ATOM%7d  %-3s %-3s  %4d    %8.3f%8.3f%8.3f%5d%8.3f 0%5.2f" % line
  rescoor = {}
  
for l in open(pdb):
  if not l.startswith("ATOM"): continue
  cres = l[21:26]
  if cres != res:
    print_res()
    res = cres
    resname = l[17:20].strip()
    if resname in mapnuc:
      if args.dna: 
        resname = mapnuc[resname][0]
      elif args.rna:
        resname = mapnuc[resname][1]
      else:
        raise ValueError("PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option" % resname)
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