# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import os

class Filename(String):
  """Data type that checks that its value is a valid local file name"""
  @staticmethod
  def typename(): return "Filename"
  def validate(self):
    from spyder.loader import file_exists
    if not file_exists(self): raise spyder.core.ValidationError
  def __new__(self, s0):
    s = s0
    if hasattr(s0, "typename") and s0.typename() == "File":
      s = s0.name
      
    l = None
    while len(s) != l:
      s = s.replace('\\\\', '\\')
      l = len(s)
    ret = String.__new__(self,s)
    Filename.validate(ret)
    return ret

spyder.__types__["Filename"] = Filename
spyder.core.arrayfunc("Filename", globals())
