# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def macro_enum(name,content):
  
  if not name.startswith("Enum"): return
  brackdepth = 0
  start = -1
  end = -1
  for n in range(len(content)):
    c = content[n]
    if c.isalnum() or c == "_": continue
    if c == "(":
      if brackdepth == 0: start = n
      brackdepth += 1
    elif c == ")":
      if brackdepth == 0:
        raise Exception("Compile error: invalid Enum member statement %s" % (name + " " + content))
      brackdepth -= 1
      if brackdepth == 0:
        end = n
        break
    elif brackdepth == 0:
      raise Exception("Compile error: invalid Enum member statement %s" % (name + " " + content))
  if start == -1 or end == -1:
    raise Exception("Compile error: invalid Enum member statement %s" % (name + " " + content))
  name, enums = content[:start], content[start+1:end]
  if enums.find(",") == -1: enums += ","
  arglist = "(" + enums + ")"
  ret = "String " + name + content[end+1:]
  ret += "\nvalidate {\n  if %s is not None: assert %s in %s\n}\nform {\n  %s.options = %s\n}\n" % (name, name, arglist,name, arglist[1:-1])
  #print ret
  return ret
  
spyder.definemacro(macro_enum)
