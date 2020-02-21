# Copyright 2007, Sjoerd de Vries
# This file is part of the Spyder module: "tarantula" 
# For licensing information, see LICENSE.txt 

class Material(String):
  """This class provides access to object materials
  See NewMaterial for more details
  """
  @staticmethod
  def typename(): return "Material"
  def __init__(self, a):
    try:
      self.construct(a)
    except:
      if len(a) > 1: raise
      self.__parse__(a[0])
  def cast(self,othertype):
    if type(othertype) == type(int): return othertype(self)
    return globals()[othertype](self)
  def __getattr__(self, method):
    if method in self.values().__dict__:
      return self.values().__dict__[method]
    else:
      m = spyder.core.method(Material, method, self, globals())
      return m
  def convert(self, target):
    c = spyder.core.convert(Material, target, self, globals())
    return c
  def construct(self, value):  
    try:
      v = get_material(value).name
    except:
      try:
        v = NewMaterial(value).name
      except:
        if not isinstance(value, str) and not isinstance(value, unicode): raise
        v = value
    String.__new__(type(self), v)
  def values(self):
    try:
      return get_material(self)
    except:
      return None
  def __parse__(self, s):
    a = spyder.core.parse(s, "Spyder")
    if len(a) == 1:
      self.construct(a)
    elif len(a) == 2:
      self.construct(a[0])


spyder.__types__["Material"] = Material
spyder.core.arrayfunc("Material", globals())
