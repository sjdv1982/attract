# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def generate_copy(typename, parentnames, source, members, deleted_members, block):
  requiredmembers, defaultmembers, optionalmembers, args, allargs = spyder.core.parse_members(typename,members,None, spyder.safe_eval)
  has_only_default = (len(requiredmembers) == 0)
  s = """  def __spydercopy__(self,_a):
    \"\"\"Auto-generated by Spyder:
     module core
     file copy.py
     function generate_copy
    Private copy constructor, for internal use only\"\"\"
"""
  if has_only_default == True:
    s += "    ok = (_a == None or isinstance(_a,%s))\n" % typename 
    #if all args are optional, at least one must be found for a non-empty arg list that is not an instance of this type
  else:
    s += "    _missing_members = []\n"
  for m in requiredmembers:
    n = m[1]
    s += "    try:\n"
    s += "      %s = _a.%s\n" % (n, n)
    s += "    except AttributeError:\n"
    s += "      _missing_members.append('%s')\n" % n
  for m in defaultmembers:
    n = m[1]
    ini = m[2]
    if m[1] in optionalmembers: ini = None
    s += "    %s = %s\n" % (n, ini)
    if has_only_default == True:
      s += "    if hasattr(_a,\"%s\") and _a.%s != None:\n      ok = True\n      %s = _a.%s\n" % (n, n, n, n)
    else: s += "    if hasattr(_a,\"%s\") and _a.%s != None: %s = _a.%s\n" % (n, n, n, n)        
  if has_only_default == True:
    s += "    if ok == False: raise spyder.ConstructionError(\"Object '%s' has no matching attributes\" % type(_a).__name__ )\n"
  else:    
    s += "    if len(_missing_members): raise spyder.ConstructionError(\"Object '%s' has missing attributes: %s\" % (type(_a).__name__, _missing_members) )\n"
  s += "    if hasattr(_a, \"__conversionstack__\"): self.__conversionstack__ = _a.__conversionstack__\n"
  s += "    return (%s,)\n" % allargs
  arraycode = """  def __spydercopy__(self,a):
    list.__init__(self, a)
    if hasattr(a, "__conversionstack__"): self.__conversionstack__ = a.__conversionstack__
 """
  if len(members) == 0: return None,arraycode  
  return s, arraycode

def generate_copydict(typename, parentnames, source, members, deleted_members, block):
  requiredmembers, defaultmembers, optionalmembers, args, allargs = spyder.core.parse_members(typename,members,None, spyder.safe_eval)
  has_only_default = (len(requiredmembers) == 0)
  s = """  def __copydict__(self,_a):
    \"\"\"Auto-generated by Spyder:
     module core
     file copy.py
     function generate_copydict
    Private dict constructor, for internal use only\"\"\"
"""
  if has_only_default == True:
    s += "    ok = (_a == None or isinstance(_a, dict))\n" #if all args are optional, at least one must be found for a non-empty arg list; empty dicts are also fine
  else:
    s += "    _missing_members = []\n"      
  for m in requiredmembers:
    n = m[1]
    s += "    try:\n"
    s += "      %s = _a['%s']\n" % (n, n)
    s += "    except KeyError:\n"
    s += "      _missing_members.append('%s')\n" % n
    s += "    except TypeError:\n"
    s += "      raise spyder.ConstructionError(\"Object '%s' does not have keys\" % type(_a).__name__)\n"      
  for m in defaultmembers:
    n = m[1]
    ini = m[2]
    if m[1] in optionalmembers: ini = None
    s += "    %s = %s\n" % (n, ini)
    if has_only_default == True:
      s += "    if \"%s\" in _a:\n      ok = True\n      %s = _a[\"%s\"]\n" % (n, n, n)
    else: s += "    if \"%s\" in _a: %s = _a[\"%s\"]\n" % (n, n, n)
    s += "    if isinstance(_a, dict) and \"%s\" in _a: %s = _a[\"%s\"]\n" % (n, n, n)                
  if has_only_default == True:
      s += "    if ok == False: raise AttributeError\n"        
  else:    
    s += "    if len(_missing_members): raise spyder.ConstructionError(\"Object '%s' has missing keys: %s\" % (type(_a).__name__, _missing_members) )\n"      
  s += "    return (%s,)\n" % allargs
  if len(members) == 0: return None,None
  return s, None

spyder.defineunimethod("__spydercopy__", generate_copy)
spyder.defineunimethod("__copydict__", generate_copydict)

