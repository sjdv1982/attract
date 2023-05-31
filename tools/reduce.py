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
  "A3": ["DA", "RA"],
  "A5": ["DA", "RA"],
  "DA3": ["DA", None],
  "DA5": ["DA", None],  
  "ADE": ["DA", "RA"],
  "C": ["DC", "RC"],
  "C3": ["DC", "RC"],
  "C5": ["DC", "RC"],
  "DC3": ["DC", None],
  "DC5": ["DC", None],
  "CYT": ["DC", "RC"],
  "G": ["DG", "RG"],
  "G3": ["DG", "RG"],
  "G5": ["DG", "RG"],
  "DG3": ["DG", None],
  "DG5": ["DG", None],  
  "GUA": ["DG", "RG"],
  "T": ["DT", None],
  "T3": ["DT", None],
  "T5": ["DT", None],
  "DT3": ["DT", None],
  "DT5": ["DT", None],  
  "THY": ["DT", None],
  "U": [None, "RU"],
  "U3": [None, "RU"],
  "U5": [None, "RU"],
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
parser.add_argument("--chain", help="Set the chain in the output PDB", default=" ")
parser.add_argument("--startres", help="Set residue number of the first residue, default 1", type=int,default=1)
parser.add_argument("--startatom", help="Set atom number of the first atom, default 1", type=int,default=1)
parser.add_argument("--batch", help="run reduce in batch mode. Input and output must be two lists of PDBs", action="store_true")


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

if args.batch and args.output is None: 
  raise ValueError("--batch requires a file list as output argument")

reducedat = os.path.split(os.path.abspath(__file__))[0] + os.sep + "reduce.dat"

def read_filelist(filelist):
  ret = []
  for l in open(filelist):
    l = l.strip()
    if not len(l): continue
    assert len(l.split()) == 1, (filelist, l)
    ret.append(l)
  return ret

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
      

def print_res(mapping=True):
  global rescounter, atomcounter, rescoor
  if not len(rescoor): return  
  rescounter += 1
  if mapping: print(res[1:].strip(), rescounter)
  for l in ff[resname]:
    if (l[0], l[1]) not in rescoor: continue
    c = rescoor[(l[0], l[1])]
    x, y, z = c[1]/c[0], c[2]/c[0], c[3]/c[0]
    if args.compat and resname in ("GLN", "GLU") and l[1] in ("CN1", "CO1"):
      x, y, z = c[1] * 0.333, c[2] * 0.333, c[3] * 0.333
    atomcounter += 1
    atomname = l[1]
    if len(atomname) < 4:
        atomname = " " + atomname + "   "[len(atomname):]
    line = (atomcounter, atomname, resname, chain, rescounter, x, y, z, l[0], l[3], 1.0)
    print("ATOM%7d %4s %3s %s%4d    %8.3f%8.3f%8.3f%5d%8.3f 0%5.2f" % line, file=outp)
  rescoor = {}
  
ff = read_forcefield(reducedat)
chain = args.chain
assert len(chain) == 1, chain

def run(pdb,mapping=True):
  global res, resname, rescounter, atomcounter, rescoor  
  res = None
  resname = None
  rescounter = args.startres-1
  atomcounter = args.startatom-1
  rescoor = {}

  
  for l in open(pdb):
    if not l.startswith("ATOM"): continue
    cres = l[21:26]
    if cres != res:
      print_res(mapping)
      res = cres
      resname = l[17:20].strip()
      if resname in mapnuc:
        if args.dna: 
          resname = mapnuc[resname][0]
        elif args.rna:
          resname = mapnuc[resname][1]
        else:
          raise ValueError("PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option" % resname)
        if resname is None:
          if args.dna: na = "DNA"
          if args.rna: na = "RNA"
          raise ValueError("'%s' can't be %s" % (l[17:20].strip(), na))      
        
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

if args.batch:
  infiles = read_filelist(args.pdb)
  for f in infiles:
    assert os.path.exists(f), f
  outfiles = read_filelist(args.output)
  for pdb, outfile in zip(infiles, outfiles):
    outp = open(outfile, "w")    
    run(pdb,mapping=False)
    print_res(mapping=False)
    outp.close()
else:  
  pdb = args.pdb
  assert os.path.exists(pdb)
  if args.output is None:
    args.output = os.path.splitext(pdb)[0] + "r.pdb"
  outp = open(args.output, "w")    
  run(pdb)
  print_res()    
