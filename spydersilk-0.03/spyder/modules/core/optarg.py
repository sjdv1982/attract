# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def macro_optarg(name, content):
  if name[0] != "*": return
  content0 = content
  for n in range(len(content)):
    if content[n].isalnum() == False and content[n] != "_":
      content = content[:n]
      break
  ret = name[1:] + " " + content0 + "\n__init__ {\n  " + content + "\n}\n"
  return ret

spyder.definemacro(macro_optarg)

