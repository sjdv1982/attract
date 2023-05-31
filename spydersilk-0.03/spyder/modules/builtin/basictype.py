# Copyright 2007, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import ast

class StringLike(Object): pass

class Float(float,Object):
  """Wrapper class around a Python float
  The conversion engine does not operate on Float"""
  @staticmethod
  def typename(): return "Float"
  def dict(self): 
    """Called by the dict function of Spyder classes
     that have Float members
    For internal use only"""
    return self
  def __eq__(self, other):
    return float(self) == other
  def __ne__(self,other): return not self.__eq__(other)
  
  def __hash__(self): return float.__hash__(self)
  def __print__(self, spaces, mode):
    """Called by the pretty-print function of Spyder classes
    that have Float members
    For internal use only"""
    return str(self)

class Integer(int,Object):
  """Wrapper class around a Python int
  The conversion engine does not operate on Integer"""  
  @staticmethod
  def typename(): return "Integer"
  def dict(self): 
    """Called by the dict function of Spyder classes
     that have Integer members
    For internal use only"""
    return self
  def __eq__(self, other):
    return int(self) == other
  def __ne__(self,other): return not self.__eq__(other)
  
  def __hash__(self): return int.__hash__(self)    
  def __print__(self, spaces, mode):
    """Called by the pretty-print function of Spyder classes
    that have Integer members
    For internal use only"""    
    return str(self)

class String(str,StringLike):
  """Wrapper class around a Python str
  The conversion engine does not operate on String, conversion into something else is always a cast"""    
  @staticmethod
  def typename(): return "String"
  def __new__(self, s):
    if s is None: raise ValueError
    if isinstance(s, String): return str.__new__(self,s) 
    s = str(s)
    if len(s) and s[0] == s[-1]:
      if s[0] in ("'", '"'):
        try:
          astree = ast.parse(s)
          s = list(ast.iter_fields(astree))[0][1][0].value.s
        except:
          pass
    ret = str.__new__(self,s)
    ret.__validate__()
    return ret
  def __validate__(self):
    pass
  def dict(self): 
    """Called by the dict function of Spyder classes
     that have String members
    For internal use only"""
    return self
  def __eq__(self, other):
    return str(self) == other
  def __ne__(self,other): return not self.__eq__(other)
 
  def __hash__(self): return str.__hash__(self)
  def __print__(self, spaces, mode):
    """Called by the pretty-print function of Spyder classes
    that have String members
    For internal use only"""        
    return str.__repr__(self)
  def convert(self, targettype):
    tar = targettype
    if isinstance(tar, str):
      tar = spyder.__types__[tar]
    return tar(self)
  def tofile(self, *args, **kargs):
    return self.convert(Data).tofile(*args, **kargs)
  def totempfile(self, *args, **kargs):
    return self.convert(Data).totempfile(*args, **kargs)
  @classmethod
  def fromfile(c, *args, **kargs):
    return c(Data.fromfile(*args, **kargs).data())

class Bool(int,Object):
  """Class that emulates a Python bool
  Unlike bool,
   "True" is equivalent to True  
   and "False" is equivalent to False
  The conversion engine does not operate on Bool"""    
  @staticmethod
  def typename(): return "Bool"  
  def __new__(self, b):
    if b == "True" or b ==  "\'True\'" or b == "\"True\"": return int.__new__(self,True)
    elif b == "False" or b == "\'False\'" or b == "\"False\"": return int.__new__(self,False)
    else: return int.__new__(self,bool(b))
  def __str__(self):
    if self == False: return "False"
    else: return "True"        
  def __repr__(self):
    if self == False: return "False"
    else: return "True"
  def dict(self): 
    """Called by the dict function of Spyder classes
     that have Bool members
    For internal use only"""
    if self: return True
    else: return False
  def __eq__(self, other):
    return bool(self) == other
  def __ne__(self,other): return not self.__eq__(other)
  
  def __hash__(self): return bool.__hash__(self)
  def __print__(self, spaces, mode):
    """Called by the pretty-print function of Spyder classes
    that have Bool members
    For internal use only"""        
    return str(self)

spyder.__types__["StringLike"] = StringLike
spyder.__types__["Integer"] = Integer
spyder.__types__["Float"] = Float
spyder.__types__["Bool"] = Bool
spyder.__types__["String"] = String

for t in ("Float", "Integer", "String", "Bool"):
  spyder.core.arrayfunc(t, globals())
