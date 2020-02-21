# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

parsed_members = {}

def parse_members(typename, members, block, safe_eval):
  if block == None and (typename,tuple(members)) in parsed_members: return parsed_members[(typename, tuple(members))]
  
  requiredmembers = [] 
  defaultmembers = []
  memberlist = set()

  optionalmembers = []
  if block != None:
    optionalmembers0 = block.split("\n")
    for o in optionalmembers0:
      o = o.strip()
      if len(o): optionalmembers.append(o)

  counter = 0
  for m in members:
    counter += 1
    m = list(m)
    spl = -1
    for n in range(0,len(m[1])):
      c = m[1][n]
      if c == "=":
        spl = n
        break
      elif c == "'" or c == "\"": break
    mm = m[:2]
    if spl > -1:
      mm = [m[0],m[1][:spl].rstrip(), m[1][spl+1:].lstrip()]
    if mm[1] in optionalmembers:
      if len(mm) == 2: mm.append("None")
    if len(mm) == 3:
      if len(mm[2].rstrip()) == 0: mm[2] = "None"
      safe_eval("x=" + mm[2])
      mm.append(counter)
      defaultmembers.append(mm)
    else:
      mm.append(counter)
      requiredmembers.append(mm)
    if mm[1] in memberlist:
      raise Exception("compile error in %s: duplicate member name %s" % (typename, mm[1]))
    memberlist.add(mm[1])
  args = ""
  allargs = ""
  for m in requiredmembers:
    args += m[1] + ","    
    allargs += m[1] + ","
  for m in defaultmembers:
    args += m[1] + "=" + m[2] + ","
    allargs += m[1] + ","
  args = args[:-1]    
  allargs = allargs[:-1]
  if block != None: parsed_members[(typename,tuple(members))] = (requiredmembers, defaultmembers, optionalmembers, args, allargs)
  return requiredmembers, defaultmembers, optionalmembers,args, allargs
