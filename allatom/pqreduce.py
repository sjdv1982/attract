import pdb2pqr
import sys, os, tempfile
import StringIO

oldsyspath = sys.path
currdir = os.path.split(__file__)[0]
if not len(currdir): currdir = "."
sys.path = [currdir]
import reduce
sys.path = oldsyspath

pdb = sys.argv[1]
transfile = sys.argv[2]
topfile = sys.argv[3]

pqrhandle, pqrfile = tempfile.mkstemp()

mapf = StringIO.StringIO()   
pdbf = StringIO.StringIO()   
pdblines = open(pdb).readlines()
reduce.run(pdblines, topfile, transfile, pdbf, mapf, [])

tmphandle, tmpfile = tempfile.mkstemp()
tmpf = open(tmpfile, "w")
for l in pdbf.getvalue().split("\n"):
  if l.find("XXX") == -1:
    print >> tmpf, l
tmpf.close()

args = [pdb2pqr.__file__, "--ff=charmm", tmpfile, pqrfile]
#if pdb2pqr.PACKAGE_PATH != "":
#  sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
oldstdout = sys.stdout
sys.stdout = sys.stderr
pdb2pqr.mainCommand(args)
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

patches = {}  
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

reduce.run(pdblines, topfile, transfile, outf, mapf2, patches)
outlines = outf.getvalue().split("\n")
for l in outlines: 
  if l.find("XXX") > -1:
    ok = False
    if l[12:16] == " HG " and  l[17:20] == "CYS": ok = True
    if l[12:16] == " OXT": ok = True
    if l[12:16] == " HT1": ok = True
    if l[12:16] == " HT2": ok = True
    if l[12:16] == " HT3": ok = True
    if not ok:
      raise ValueError("Missing atom:\n%s" % l)

outfile = os.path.splitext(pdb)[0] + "-aa.pdb"
out = open(outfile, "w")
for l in outlines:
  if len(l): print >> out, l

print mapf.getvalue()

