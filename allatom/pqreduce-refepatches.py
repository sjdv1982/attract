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
pdbref = sys.argv[4]

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
  #sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
oldstdout = sys.stdout
sys.stdout = sys.stderr
pdb2pqr.mainCommand(args)
sys.stdout = oldstdout
pqr = os.fdopen(pqrhandle)
pqrlines = pqr.readlines()
pqr.close()
os.remove(pqrfile)

pqrrefhandle, pqrreffile = tempfile.mkstemp()
mapfref = StringIO.StringIO()   
pdbfref = StringIO.StringIO()   
pdbreflines = open(pdbref).readlines()
reduce.run(pdbreflines, topfile, transfile, pdbfref, mapfref, [])
tmphandle, tmpfile = tempfile.mkstemp()
tmpf = open(tmpfile, "w")
for l in pdbfref.getvalue().split("\n"):
  if l.find("XXX") == -1:
    print >> tmpf, l
tmpf.close()

args = [pdb2pqr.__file__, "--ff=charmm", tmpfile, pqrreffile]
#if pdb2pqr.PACKAGE_PATH != "":
  #sys.path.extend(pdb2pqr.PACKAGE_PATH.split(":"))
oldstdout = sys.stdout
sys.stdout = sys.stderr
pdb2pqr.mainCommand(args)
sys.stdout = oldstdout
pqrref = os.fdopen(pqrrefhandle)
pqrreflines = pqrref.readlines()
pqrref.close()
os.remove(pqrreffile)

pdblines = []
pdbreflines = []
his = {}

repl = (
  (" H  ", " HN "),
  (" H  ", " HT1"),
  (" H2 ", " HT2"),
  (" H3 ", " HT3"),
)  
  
counter = max(len(pqrreflines),len(pqrlines))
for lnr in range(counter):
  if lnr < len(pqrlines):
    l = pqrlines[lnr]
    pdblines.append(l)
    atom2 = l[12:16]
    for pin, pout in repl:
      if atom2 == pin:
	p = l[:12] + pout + l[16:]
	pdblines.append(p)
      
  if lnr < len(pqrreflines):
    lref = pqrreflines[lnr]
    pdbreflines.append(lref)
    atom2 = lref[12:16]
    for pin, pout in repl:
      if atom2 == pin:
	p = lref[:12] + pout + lref[16:]
	pdbreflines.append(p)
      
    atom = lref[13:16]
    resname = lref[17:20]
    resid = lref[21:27]  
    if resname == "HIS":
      if resid not in his: his[resid] = []
      if atom in ("HE1", "HE2", "HD1", "HD2"): 
	his[resid].append(atom)
	print resid,atom
      

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
    if l[12:16] == " HT3" and  l[17:20] == "PRO": ok = True
    if 'HIS' in l: ok=True
    if not ok:
      raise ValueError("Missing atom:\n%s" % l)

outfile = os.path.splitext(pdb)[0] + "-aa.pdb"
out = open(outfile, "w")
for l in outlines:
  if len(l): print >> out, l
  
out.close()
  
outfref = StringIO.StringIO()   
mapf2ref = StringIO.StringIO()   

reduce.run(pdbreflines, topfile, transfile, outfref, mapf2ref, patches)
outlines = outfref.getvalue().split("\n")
for l in outlines: 
  if l.find("XXX") > -1:
    ok = False
    if l[12:16] == " HG " and  l[17:20] == "CYS": ok = True
    if l[12:16] == " HT3" and  l[17:20] == "PRO": ok = True
    if not ok:
      raise ValueError("Missing atom:\n%s" % l)

outfile = os.path.splitext(pdbref)[0] + "-aa.pdb"
out = open(outfile, "w")
for l in outlines:
  if len(l): print >> out, l

out.close()
print mapf.getvalue()

