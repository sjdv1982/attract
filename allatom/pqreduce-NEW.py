"""
Reduce converts a PDB into a modified PDB format, where the atom type for each atom is indicated
This reduced PDB is what is understood by ATTRACT.
These atoms are derived from the transfile and topology file 
"""

import pdb2pqr
import sys, os, tempfile
import StringIO

has_argparse = False
try:
  import argparse  
  has_argparse = True  
except ImportError:
  import optparse  #Python 2.6

oldsyspath = sys.path
currdir = os.path.split(__file__)[0]
if not len(currdir): currdir = "."
sys.path = [currdir]
import reduce 
sys.path = oldsyspath

if has_argparse:
  parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("pdb",help="PDB file to reduce")
  parser.add_argument("transfile",help="Trans file that contains the atom types")
  parser.add_argument("topfile",help="Topology file in CNS format")
  parser.add_argument("output",help="all-atom reduced output PDB file", nargs="?")
else:
  parser = optparse.OptionParser()
  parser.add_argument = parser.add_option
parser.add_argument("--notermini",help="Do not add N- and C-termini for each chain", action="store_true")
parser.add_argument("--refe",help="Reference file to determine histidine/cysteine states")
parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
parser.add_argument("--whatif",help="Uses the WHATIF server instead of pqreduce (automatically enabled with --rna or --dna)", action="store_true")

if has_argparse:
  args = parser.parse_args()
else:
  args, positional_args = parser.parse_args()
  args.pdb = None
  args.output = None
  if positional_args:
    args.pdb = positional_args[0]
    args.transfile = positional_args[1]
    args.topfile = positional_args[2]
    if len(positional_args) > 3: args.output = positional_args[3]

pqrhandle, pqrfile = tempfile.mkstemp()

mapf = StringIO.StringIO()   
pdbf = StringIO.StringIO()   
pdblines = open(args.pdb).readlines()
reduce.run(pdblines, args.topfile, args.transfile, pdbf, mapf, [])

tmphandle, tmpfile = tempfile.mkstemp()
tmpf = open(tmpfile, "w")
for l in pdbf.getvalue().split("\n"):
  if l.find("XXX") == -1:
    print >> tmpf, l
tmpf.close()

patches = {}  

if args.whatif or args.rna or args.dna:
  pqrlines = [] ###TODO
else:
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
  os.remove(pqrfile)

pdblines = []
his = {}

repl = (
  (" H  ", " HN "),
  (" H  ", " HT1"),
  (" H2 ", " HT2"),
  (" H3 ", " HT3"),
)  
for lnr in range(len(pqrlines)):
  l = pqrlines[lnr]
  pdblines.append(l)
  atom2 = l[12:16]
  for pin, pout in repl:
    if atom2 == pin:
      p = l[:12] + pout + l[16:]
      pdblines.append(p)
  atom = l[13:16]
  resname = l[17:20]
  resid = l[21:27]  
  if resname == "HIS":
    if resid not in his: his[resid] = []
    if atom in ("HE1", "HE2", "HD1", "HD2"): his[resid].append(atom)

for resid in his:
  protons = his[resid]
  if len(protons) == 4: continue
  assert len(protons) == 3, (resid, protons)
  if "HE2" not in protons or "HE1" not in protons: 
    patches[resid] = "hisd"
  elif "HD2" not in protons or "HD1" not in protons: 
    patches[resid] = "hise"
  else:  
    raise ValueError((resid, protons))

outf = StringIO.StringIO()   
mapf2 = StringIO.StringIO()   

termini = True
if args.notermini: termini = False
reduce.run(pdblines, topfile, transfile, outf, mapf2, patches, termini=termini)
outlines = outf.getvalue().split("\n")
for l in outlines: 
  if l.find("XXX") > -1:
    ok = False
    if l[12:16] == " HG " and  l[17:20] == "CYS": ok = True
    if termini:
      if l[12:16] == " OXT": ok = True
      if l[12:16] == " HT1": ok = True
      if l[12:16] == " HT2": ok = True
      if l[12:16] == " HT3": ok = True
    if not ok:
      raise ValueError("Missing atom:\n%s" % l)

if args.output is None:
  args.output = os.path.splitext(pdb)[0] + "-aa.pdb"
out = open(args.output, "w")
for l in outlines:
  if len(l): print >> out, l

print mapf.getvalue()

