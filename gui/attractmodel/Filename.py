# Copyright 2007-2012, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 
#
# Changes to builtin/Filename.py, will be part of Spyder 3.3
#

import os

class Filename(String):
  """Data type that checks that its value is a valid local file name"""
  @staticmethod
  def typename(): return "Filename"
  def validate(self):
    if not os.path.exists(self): raise spyder.core.ValidationError
  def __new__(self, s0):
    s = s0
    if isinstance(s0, tuple) and len(s0) == 2:
      s = s0[0]
    if hasattr(s0, "typename") and s0.typename() == "File":
      s = s0.name      
    l = None
    while len(s) != l:
      s = s.replace('\\\\', '\\')
      l = len(s)
    ret = String.__new__(self,s)
    Filename.validate(ret)
    return ret
  name = property(lambda self: self)
  value = property(lambda self: self)
  
spyder.__types__["Filename"] = Filename
spyder.core.arrayfunc("Filename", globals())
