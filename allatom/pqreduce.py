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

args = [pdb2pqr.__file__, "--ff=charmm", pdb, pqrfile]
if pdb2pqr.PACKAGE_PATH != "":
  sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
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
for lnr in range(len(pqrlines)):
  l = pqrlines[lnr]
  l = l.replace(" H ", " HN")
  resid = l[21:27]
  if l.find(" H2 ") > -1:
    for ll in pqrlines[lnr+1:lnr+10]: print l,
    for n in range(len(pdblines)):
      if pdblines[n][21:27] == resid:
        pdblines[n] = pdblines[n].replace( " HN ", " HT1")
    for n in range(lnr+1,len(pqrlines)):
      if pqrlines[n][21:27] == resid:
        pqrlines[n] = pqrlines[n].replace( " H  ", " HT1")
    l = l.replace(" H2 ", " HT2")  
  l = l.replace(" H3 ", " HT3")
  pdblines.append(l)
  atom = l[13:16]
  resname = l[17:20]
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
mapf = StringIO.StringIO()   
reduce.run(pdblines, topfile, transfile, outf, mapf, patches)
outlines = outf.getvalue().split("\n")
for l in outlines: 
  if l.find("XXX") > 1:
    raise ValueError("Missing atom:\n%s" % l)

outfile = os.path.splitext(pdb)[0] + "-aa.pdb"
out = open(outfile, "w")
for l in outlines:
  if len(l): print >> out, l

print mapf.getvalue()

