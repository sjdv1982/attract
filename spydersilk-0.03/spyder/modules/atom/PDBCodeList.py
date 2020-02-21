## DEFINE PDBCodeList, PDBCodeChainList

class PDBCodeList(PDBCodeArray):
  def __parse__(self, s):
    try:
      if s.startswith('"') and s.endswith('"'): s = s[1:-1]
      if s.startswith("'") and s.endswith("'"): s = s[1:-1]    
    except:
      pass
    try:
      ret = PDBCodeArray.__parse__(self, s)
    except:      
      ret = PDBCodeArray.__parse__(self, "[" + s + "]")
    return ret
  def __spydercopy__(self, *args, **kargs):
    if len(args):
      v = None
      if isinstance(args[0],str): v = args[0]
      elif isinstance(args[0],tuple) and len(args[0]): v = args[0][0]
      if v != None and v.find(",") > -1:
        raise ValueError("Invalid PDB code %s" % v)
    return PDBCodeArray.__spydercopy__(self, *args, **kargs)
  def __repr__(self, *args, **kargs):
    return "'" + self.__str__(*args, **kargs) + "'"
  def __str__(self, *args, **kargs):
    return ",".join([str(v) for v in self])
  def __print__(self, *args, **kargs):
    return self.__repr__(*args,**kargs)

spyder.__types__["PDBCodeList"] = PDBCodeList
spyder.core.arrayfunc("PDBCodeList", globals())

class PDBCodeChainList(PDBCodeChainArray):
  def __parse__(self, s):
    try:
      if s.startswith('"') and s.endswith('"'): s = s[1:-1]
      if s.startswith("'") and s.endswith("'"): s = s[1:-1]    
    except:
      pass
    try:
      ret = PDBCodeChainArray.__parse__(self, s)
    except:      
      ret = PDBCodeChainArray.__parse__(self, "[" + s + "]")
    return ret
  def __spydercopy__(self, *args, **kargs):
    if len(args):
      v = None
      if isinstance(args[0],str): v = args[0]
      elif isinstance(args[0],tuple) and len(args[0]): v = args[0][0]
      if v != None and v.find(",") > -1:
        raise ValueError("Invalid PDB code %s" % v)
    return PDBCodeChainArray.__spydercopy__(self, *args, **kargs)    
  def __repr__(self, *args, **kargs):
    return "'" + self.__str__(*args, **kargs) + "'"
  def __str__(self, *args, **kargs):
    return ",".join([str(v) for v in self])
  def __print__(self, *args, **kargs):
    return self.__repr__(*args,**kargs)

spyder.__types__["PDBCodeChainList"] = PDBCodeList
spyder.core.arrayfunc("PDBCodeChainList", globals())
