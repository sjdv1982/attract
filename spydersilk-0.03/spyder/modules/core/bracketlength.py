# Copyright 2007-2013, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def macro_bracketlength(name,content): 
  content0 = content
  default = ""
  p = content0.find("=")
  if p > -1:
    content0 = content[:p].rstrip()
    default = content[p:]
  if not content0.endswith("]"): return
  p = content0.rfind("[")
  ret = name + " " + content0[:p] + default
  l = content0[p+1:-1]
  v = content0[:p]
  ret += "\nvalidate {\n  assert %s == None or len(%s) == %s\n}\nform {\n  %s.length = %s\n  %s.form = \"hard\"\n}\n" % (v,v, l,v,l,v)
  return ret

spyder.definemacro(macro_bracketlength)
