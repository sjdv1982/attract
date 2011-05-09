def read_struc(fil):
  lines = open(fil).readlines()
  header = []
  mode = 0
  structures = []
  centeredlen = 0
  for l in lines:
    l = l.rstrip("\n")
    if centeredlen < 2:
      if l.startswith("##"): continue
      header.append(l)
      if l.startswith("#centered"): centeredlen += 1
      continue
    if mode == 2 and l[0] == "#": mode = 0
    if mode == 0:
      assert int(l[1:]) == len(structures)+1
      mode = 1
      structures.append(([],[]))
      continue
    if mode == 1:
      if l[:2] == "##":
	structures[-1][0].append(l)
      else:
	mode = 2
    if mode == 2:
      structures[-1][1].append(l)
  return header, structures
