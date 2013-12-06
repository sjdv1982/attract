import sys, os, glob

runname = sys.argv[1]
peptidebond = float(sys.argv[2])   #should be about 1.4 A, but can be increased to reflect uncertainty
residuelength = float(sys.argv[3]) #should be about 3.4-3.8 A
patterns = sys.argv[4:]

assert len(patterns) > 0

def get_resids(pdbfile):
  r = None
  resids = []
  for l in open(pdbfile).readlines():
    if not l.startswith("ATOM"): continue
    resid = l[21:26]
    if resid != r:
      r = resid
      resids.append(r)
  return resids

def get_links(linkfile):
  links = []
  lines = open(linkfile).readlines()
  for l in lines:
    link = tuple([int(v) for v in l.split()])
    links.append(link)
  return links


def find_atom(data, resnr, atom): 
  for lnr,l in enumerate(data):
    #print(l)
    if l[:4] != "ATOM": continue
    r = int(l[22:26])
    a = l[13:16]
    #print(r, a)
    if r == resnr and a == atom: return lnr
  raise ValueError((resnr,atom))
  
def process_pattern(chain, pattern, pdblines, offset, restsele, restdef):
  master_pdb = pattern + ".pdb"
  resids = get_resids(master_pdb)

  links = get_links(pattern + ".links")
  
  bpdbs = glob.glob(pattern + ".pdb-[0-9]") + glob.glob(pattern + ".pdb-?[0-9]") \
   + glob.glob(pattern + ".pdb-??[0-9]")
  
  reduced_data = []
  bresids = []
  for bpdb in bpdbs:
    bresids0 = get_resids(bpdb)
    bresids0 = [resids.index(l) for l in bresids0]
    bresids.append(bresids0)
    
    cmd = "\cp %s /tmp/ssetup.pdb" % bpdb
    os.system(cmd)
    cmd = "$ATTRACTDIR/reduce /tmp/ssetup.pdb > /dev/null; \cp /tmp/ssetupr.pdb %sr" % bpdb
    os.system(cmd)
    rdata = open("%sr" % bpdb).readlines()    
    reduced_data.append(rdata)

  boffsets = []
  for r in reduced_data:
    pdblines += r
    boffsets.append(offset)
    offset += len(r)
    pdblines.append("TER\n")

  for body1, body2, resid1, resid2, linkdist in links:
    id1 = bresids[body1-1]
    id2 = bresids[body2-1]
    d1 = reduced_data[body1-1]
    d2 = reduced_data[body2-1]
    resnr1 = id1.index(resid1-1) + 1
    resnr2 = id2.index(resid2-1) + 1
    offset1 = boffsets[body1-1]
    offset2 = boffsets[body2-1]
    atom1 = find_atom(d1, resnr1, "C  ") + offset1 + 1
    atom2 = find_atom(d2, resnr2, "N  ") + offset2 + 1
    print atom1, atom2
    sele1 = "%s_body%d_%d_C" % (chain, body1, int(resids[resid1-1][1:]))
    sele2 = "%s_body%d_%d_N" % (chain, body2, int(resids[resid2-1][1:]))
    restsele.append("%s 1 %d" % (sele1, atom1))
    restsele.append("%s 1 %d" % (sele2, atom2))
    dist = peptidebond + residuelength * linkdist
    restdef.append("%s %s 1 %.3f 1000" % (sele1, sele2, dist))
    
    
  return offset, len(bpdbs)
  
offset = 0
bodies = 0
pdblines = []
restsele = [] #restraint selections
restdef = [] #restraint definitions
for pnr,pat in enumerate(patterns):
  ch = chr(64+pnr+1)
  offset, pbodies = process_pattern(ch, pat, pdblines, offset, restsele, restdef)
  bodies += pbodies

f = open(runname +"r.pdb", "w")
f.write("%d\n" % bodies)
for l in pdblines:
  f.write(l)
f.write("END\n")
f.close()

f = open(runname +".rest", "w")
for l in restsele:
  f.write(l + "\n")
f.write("\n")  
for l in restdef:
  f.write(l + "\n")
f.close()
