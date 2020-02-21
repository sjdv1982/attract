# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import functools
import Spyder

class ObjectList(spyder.core.spyderlist):
  """This class is similar to a Spyder Array, but it its members need not
  be of a single Spyder class.
  Instead, every member must be of any Spyder class, i.e. it must be an instance of Spyder.Object
  """
  __constructor__ = "constructor_fromany"
  
  @staticmethod
  def typename(): return "ObjectList"
  def __init__(self, *args, **kargs):
    return ObjectList.__dict__[ObjectList.__constructor__](self, *args, **kargs)
  def constructor_fromany(self, *a):
    try:
      spyder.core.spyderlist.__init__(self,*a)
      self.validate()
    except Exception as e:
      if len(a) > 1: raise
      self.__parse__(a[0])
      self.validate()
  def constructor_fromdict(self, dic):      
    if isinstance(dic, dict):
      aa = []
      for key in dic:
        try:
          constr = globals()[key].fromdict          
        except AttributeError:
          constr = globals()[key]
        aa.append(constr(dic[key]))
      spyder.core.spyderlist.__init__(self,*aa)
      self.validate()    
    elif isinstance(dic, list):      
      spyder.core.spyderlist.__init__(self,*dic)
      self.validate()          
  def cast(self,othertype):
    if type(othertype) == type(int): return othertype(self)
    return globals()[othertype](self)
  def __getattr__(self, method):
    m = spyder.core.method(ObjectList, method, self)
    return m
  def convert(self, target):
    c = spyder.core.convert(ObjectList, target, self)
    return c
  def validate(self):
    try:
      self.__validate__()
    except Exception as inst: 
      if isinstance(inst, spyder.core.ValidationError):
        raise 
      else: raise spyder.core.ValidationError(inst)  
  def __validate__(self):
    for o in self:
      assert isinstance(o, Spyder.Object), (type(o), o)
  def list(self):
    return list(self)
  def length(self): return len(self)
  def __parse__(self, s):
    a = spyder.core.parse(s, "Spyder")
    if len(a) == 1:
      spyder.core.spyderlist.__init__(self,a)
    elif len(a) == 2:
      spyder.core.spyderlist.__init__(self,a[0])
  def __print__(self,spaces,mode):
    ret = "%s (\n" % self.typename()
    for o in self:
      p = o.__print__(spaces+2,mode)
      tp = o.typename()
      if tp in ("Integer","Float","Bool","String"): p = tp+"("+p+")"
      ret += (spaces+2) * " " + p + ",\n"
    ret += spaces * " " + ")"
    return ret
  def __str__(self): return self.__print__(0, "str")
  def str(self): return self.__print__(0, "str")
  def data(self): return str(self)
  def __repr__(self): return self.__print__(0, "repr")
  def repr(self): return self.__print__(0, "repr")
  @classmethod
  def fromfile(c, filename, fastparse=False):
    raise NotImplementedError
  def tofile(self, filename):
    raise NotImplementedError
  def totempfile(self):
    raise NotImplementedError    
  def dict(self): 
    """Called by the dict function of Spyder classes
     that have ObjectList members
    For internal use only"""
    return self

def objectlistshow(ol,*args,**kargs):
  oldstack = spyder.core.__conversionstack__.l  
  spyder.core.conversionstack.l = []
  for o in ol:
    o.show(*args,**kargs)
  spyder.core.conversionstack.l = oldstack
  
spyder.core.definemethod("show","ObjectList", objectlistshow)
del objectlistshow

spyder.core.defineconverter("ObjectList", "None", "SPLIT")

spyder.__types__["ObjectList"] = ObjectList
spyder.core.arrayfunc("ObjectList", globals())

ObjectList.fromdict = functools.partial(spyder.__constructor, 
  "constructor_fromdict",
  ObjectList,
  "constructor_fromany",
)
