class FloatList(FloatArray):
  def __arrayvalidate__(self):    
    for v in range(len(self)): 
      try:
        Float(self[v])
      except:
        raise ValueError("Invalid number '%s'. Please specify as comma-separated numbers, e.g. 1.8, 2.7 , 3.0" % self[v])
    FloatArray.__arrayvalidate__(self)
  def __parse__(self, s):
    try:
      if s.startswith('"') and s.endswith('"'): s = s[1:-1]
      if s.startswith("'") and s.endswith("'"): s = s[1:-1]    
    except:
      pass
    try:
      ret = FloatArray.__parse__(self, s)
    except:      
      ret = FloatArray.__parse__(self, "[" + s + "]")
    return ret
  def __repr__(self, *args, **kargs):
    return "'" + self.__str__(*args, **kargs) + "'"
  def __str__(self, *args, **kargs):
    return ", ".join([str(v) for v in self])
  def __print__(self, *args, **kargs):
    return self.__repr__(*args,**kargs)
  def __eq__(self, other):
    if isinstance(other, str):
      return other == str(self)
    return FloatArray.__eq__(self, other)

spyder.__types__["FloatList"] = FloatList
spyder.core.arrayfunc("FloatList", globals())
