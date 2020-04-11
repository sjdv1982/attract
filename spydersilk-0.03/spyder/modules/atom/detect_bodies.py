def detect_bodies(pdb):
  cutoffsq1 = 5 * 5
  cutoffsq2 = 5 * 5
  lines = [l for l in pdb.splitlines() if l.startswith("ATOM")]
  pairs = []
  currpairnr = None
  def match_pair(pairnr, begin, cutoffsq):    
    z1 = (len(pairs[pairnr][0]) == 0)
    z2 = (len(pairs[pairnr][1]) == 0)
    if z1 or z2: return False
    if begin:       
      coor = pairs[pairnr][0]
      ends = [p[1] for p in pairs]      
    else: 
      coor = pairs[pairnr][1]
      ends = [p[0] for p in pairs]            
    for enr, e in enumerate(ends):
      if enr == pairnr: continue
      zz1 = (len(pairs[enr][0]) == 0)
      zz2 = (len(pairs[enr][1]) == 0)          
      if zz1 or zz2: return False      
      dissq = (e[0]-coor[0])**2 + (e[1]-coor[1])**2 + (e[2]-coor[2])**2
      if dissq < cutoffsq:
        if begin:
          pairs[enr][1] = pairs[pairnr][1]          
        else:
          pairs[enr][0] = pairs[pairnr][0]
        pairs[enr][2] += pairs[pairnr][2]
        pairs.pop(pairnr)
        return True
    return False
  oldresnr = None
  for l in lines:
    atom = l[13:16]
    resnr = int(l[22:26])
    if resnr != oldresnr:
      if currpairnr != None:
        pairs.pop()
        currpairnr = None
      oldresnr = resnr
    if atom == "N  ":
      coor = (float(l[30:38]), float(l[38:46]), float(l[46:54]))      
      if currpairnr == None:
        currpairnr = len(pairs)
        pairs.append([coor, []])
      elif len(pairs[currpairnr]) == 2:
        pairs[currpairnr][0] = coor
        pairs[currpairnr].append([resnr])
        match_pair(currpairnr, True, cutoffsq1)
        currpairnr = None
    if atom == "C  ":
      coor = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
      if currpairnr == None:
        currpairnr = len(pairs)
        pairs.append([[],coor])
      elif len(pairs[currpairnr]) == 2:
        pairs[currpairnr][1] = coor      
        pairs[currpairnr].append([resnr])
        match_pair(currpairnr, True, cutoffsq1)
        currpairnr = None  
  found = True
  while found and len(pairs) > 1:
    found = False
    for n in range(len(pairs)):
      if match_pair(n, True, cutoffsq2):
        found = True
        break
      if match_pair(n, False, cutoffsq2):
        found = True
        break
  return [p[2] for p in pairs if len(p) == 3 and len(p[0]) == 3 and len(p[1]) == 3]
