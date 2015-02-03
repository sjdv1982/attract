class HaddockResidueList(IntegerArray):
  def __arrayvalidate__(self):    
    for v in range(len(self)): 
      try:
        Integer(self[v])
      except:
        raise ValueError("Invalid residue number '%s'. Please specify your residues as comma-separated numbers, e.g. 1,2,3." % self[v])
    IntegerArray.__arrayvalidate__(self)
  def __parse__(self, s):
    try:
      if s.startswith('"') and s.endswith('"'): s = s[1:-1]
      if s.startswith("'") and s.endswith("'"): s = s[1:-1]    
    except:
      pass
    try:
      ret = IntegerArray.__parse__(self, s)
    except:      
      ret = IntegerArray.__parse__(self, "[" + s + "]")
    return ret
  def __repr__(self, *args, **kargs):
    return "'" + self.__str__(*args, **kargs) + "'"
  def __str__(self, *args, **kargs):
    return ",".join([str(v) for v in self])
  def __print__(self, *args, **kargs):
    return self.__repr__(*args,**kargs)
  def __eq__(self, other):
    if isinstance(other, str):
      return other == str(self)
    return IntegerArray.__eq__(self, other)

spyder.__types__["HaddockResidueList"] = HaddockResidueList
spyder.core.arrayfunc("HaddockResidueList", globals())
