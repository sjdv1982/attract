#import pdb2pqr
#import vmd
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

vmdhandle, vmdfile = tempfile.mkstemp()

mapf = StringIO.StringIO()   
pdbf = StringIO.StringIO()   
vmdlines = open(pdb,'r').readlines()

#reduce.run(pdblines, topfile, transfile, pdbf, mapf, [], termini=False)
#TODO: 5PHO/5TER and 3TER patches?? DNA patches (now RNA by default) => adapt reduce.py

tmphandle, tmpfile = tempfile.mkstemp()
tmpf = open(tmpfile, "w")
for l in pdbf.getvalue().split("\n"):
  if l.find("XXX") == -1:
    print >> tmpf, l
tmpf.close()

pdblines = []
his = {}

repl = (
  ("H2''", "HO2'"),
)  
for lnr in range(len(vmdlines)):
  l = vmdlines[lnr]
  lsplit=l.split()  					# added by Isaure
  if len(lsplit)==0:continue 
  if lsplit[0]!='ATOM':continue
  if lsplit[-1]=='H' and toremove: continue  
  if lsplit[-1]!='H':
	if lsplit[-4]=='0.00': 
		toremove=True
		continue  				# added by Isaure
	else:
		toremove=False
  resname = l[17:20]
  if resname == "URA":
    l = l[:17] + "URI" + l[20:]
  pdblines.append(l)
  atom2 = l[12:16]
  for pin, pout in repl:
    if atom2 == pin:
      p = l[:12] + pout + l[16:]
      pdblines.append(p)
  atom = l[13:16]
  resid = l[21:27]

outf = StringIO.StringIO()   
mapf2 = StringIO.StringIO()   

reduce.run(pdblines, topfile, transfile, outf, mapf2, {}, termini=False)
outlines = outf.getvalue().split("\n")
reoutlines=[]
for l in outlines: 
  if l.find("XXX") > -1: continue
  reoutlines.append(l)

outfile = os.path.splitext(pdb)[0] + "-aa.pdb"
out = open(outfile, "w")
for l in reoutlines:
  if len(l): print >> out, l

print mapf.getvalue()
