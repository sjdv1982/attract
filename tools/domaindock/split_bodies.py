margin = 4 #number of residues to be known beyond SS element boundaries
minsize_assign = 3 #minimum amount of identical consecutive SS assignments for a SS element
minsize_final = 9 #minimum size of a SS element

assert 2*margin+1 >= minsize_final

resids = []
data = {}

def parse_pdb(pdbfile):
  r = None
  for l in open(pdbfile).readlines():
    if not l.startswith("ATOM"): continue
    resid = l[21:26]
    if resid != r:
      r = resid
      assert r not in data #contiguous definition of a residue
      data[r] = []
      resids.append(r)
    data[r].append(l)  

def parse_dssp(dsspfile):
  lines0 = open(dsspfile).readlines()
  lines = []
  parse = False
  for l in lines0:  
    if not len(l.strip()): continue
    if not parse and l.startswith("  #  RESIDUE AA STRUCTURE BP1 BP2"):
      parse = True
      continue
    if not parse: continue
    lines.append(l)
  if len(lines) != len(resids):
    raise Exception("Different number of residues in DSSP (%d) vs PDB (%d)" % (len(lines), len(resids)))
  
  bodymap = [None] * len(resids)

  sheet_to_bodies = {}
  curr_ss = None
  bodycount = 0
  for lnr,l in enumerate(lines):
    resid = l[11] + l[6:10]
    if resid != resids[lnr]:
      raise Exception("Different identifier for residue %d in DSSP ('%s') vs PDB ('%s')" % (resid, resids[lnr]))
    ss = l[16]
    if ss == curr_ss:
      if ss in (" ","B","T","S"): continue
      bodymap[lnr] = curr_body
      continue
      
    #new element
    curr_ss = ss
    sheet = l[33]
    if sheet != " ":
      if sheet in sheet_to_bodies: 
        curr_body = sheet_to_bodies[sheet]
      else:
        curr_body = bodycount  
        sheet_to_bodies[sheet] = curr_body
        bodycount+=1      
    else:      
      if ss in (" ","B","T","S"): continue
      curr_body = bodycount
      bodycount += 1

    bodymap[lnr] = curr_body
    
  return bodycount, bodymap
 
 
def filter_bodymap(bodycount, bodymap, minlen):
  counts = [0] * bodycount
  for l in bodymap: 
    if l is not None: counts[l] += 1
  
  newbodycount = bodycount
  for body in reversed(range(bodycount)):
    if counts[body] < minlen:
      newbodycount -= 1
      for lnr,l in enumerate(bodymap): 
        if l == body:  
          bodymap[lnr] = None
        elif l > body:  
          bodymap[lnr] -= 1
  
  return newbodycount
      
def get_inter_segments(bodymap):
  segments = []
  prev = None
  curr_none = True
  start = 0
  for lnr, l in enumerate(bodymap):
    if curr_none:
      if l is None: continue
      currlen = lnr - start
      if currlen > 0:
        segments.append((prev, l, start+1, currlen))
      prev = l
      curr_none = False
    else:
      if l is not None: 
        prev = l      
        continue
      start = lnr
      curr_none = True  
  if curr_none: segments.append((prev, None, start+1, len(bodymap)-start))
  return segments

def fill_bodymap(bodymap, maxlen):
  segments = get_inter_segments(bodymap)
  for prev, next, start, length in segments:
    if length > maxlen: continue
    if prev is None and next is None: continue
    v1 = prev if prev is not None else next
    v2 = next if next is not None else prev
    if v1 != v2: continue
    for n in range(length):
      bodymap[start-1+n] = v1

def interpolate_bodymap(bodymap, margin):    
  segments = get_inter_segments(bodymap)
  for prev, next, start, length in segments:
    if prev is None and length <= margin:
      for n in range(length):
        bodymap[start-1+n] = next
    elif next is None and length <= margin:
      for n in range(length):
        bodymap[start-1+n] = prev
    elif length <= 2*margin:
      mid = int(length/2.0+0.999) #rounded up
      for n in range(length):
        v = next if n+1 > mid else prev
        bodymap[start-1+n] = v

def fill_loops(bodycount, bodymap):
  segments = get_inter_segments(bodymap)
  for prev, next, start, length in segments:
    for n in range(length):
      bodymap[start-1+n] = bodycount
    bodycount += 1  
  return bodycount

def get_links(bodymap):
  links = []
  prev = None
  prevnr = None
  for lnr, l in enumerate(bodymap):
    if l is None: continue
    if l == prev: 
      prevnr = lnr
      continue
    if prev is not None:
      links.append((prev, l, prevnr, lnr, lnr-prevnr-1))
    prev = l
    prevnr = lnr
  return links       

def write_bodies(pdbfile, data, resids, bodycount, bodymap):
  files = []
  for n in range(bodycount):
    f = open(pdbfile + "-" + str(n+1), "w")
    files.append(f)
  for lnr, l in enumerate(bodymap):
    resid = resids[lnr]
    d = data[resid]
    f = files[l]
    for atom in d:
      f.write(atom)
  for f in files: f.close()  
  
      
import sys
pdbfile = sys.argv[1]
dsspfile = sys.argv[2]

parse_pdb(pdbfile)
bodycount, bodymap = parse_dssp(dsspfile)

bodycount = filter_bodymap(bodycount, bodymap, minsize_assign)
fill_bodymap(bodymap, 2*margin)
interpolate_bodymap(bodymap, margin)
bodycount = filter_bodymap(bodycount, bodymap, minsize_final)
fill_bodymap(bodymap, 2*margin)

###############################################################
#Now all we have left unassigned is long loops (length > 2 * margin)

#For now, treat them as separate bodies
bodycount = fill_loops(bodycount, bodymap)

#Alternatively, we could nibble off the margins and treat the rest as missing atoms
#extrapolate_bodymap(bodymap, margin)

###############################################################

links = get_links(bodymap)
for link in links: 
  print link[0]+1, link[1]+1, link[2]+1, link[3]+1, link[4]

write_bodies(pdbfile,  data, resids, bodycount,  bodymap)

