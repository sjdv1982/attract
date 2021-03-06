# Copyright 2009-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def generate_eq(typename, parentnames, source, members, deleted_members, block):
  requiredmembers, defaultmembers, optionalmembers, args, allargs = spyder.core.parse_members(typename,members,None, spyder.safe_eval)
  s = """  def __ne__(self,_a): return not self.__eq__(_a)
  def __eq__(self,_a):
    \"\"\"Auto-generated by Spyder:
     module core
     file eq.py
     function generate_eq
    Comparison operator\"\"\"
"""
  arraycode = s
  s += "    try:\n"
  s += "      if not isinstance(_a,%s) and not isinstance(self, type(_a)): return False\n" % typename
  s += "    except (TypeError, AttributeError, spyder.ValidationError): return False\n"
  s += "    if _a is self: return True\n"
  for m in requiredmembers + defaultmembers:
    s += "    if self.%s != _a.%s: return False\n" % (m[1], m[1])
  if len(parentnames):
    for parentname in parentnames:
      s += "    try:\n"
      s += "      if not %s.__eq__(self, _a): return False\n" % parentname
      s += "    except AttributeError:\n"
      s += "      pass\n"
  s += "    return True\n"
  arraycode += "    if not isinstance(_a, %(classname)s) and not isinstance(self, type(_a)): return False\n" 
  arraycode += "    if _a is self: return True\n"
  arraycode += "    if len(self) != len(_a): return False\n"
  arraycode += "    for n,nn in zip(self,_a):\n"
  arraycode += "      if n != nn: return False\n"
  arraycode += "    return True\n"
  return s, arraycode

spyder.defineunimethod("__eq__", generate_eq)
