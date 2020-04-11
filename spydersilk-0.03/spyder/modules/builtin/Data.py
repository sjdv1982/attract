# Copyright 2007,2013, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import ast, ctypes

transtable = ""
for n in range(0,256):
  if chr(n) in ('\n', '\t', '\r'):
    transtable += chr(n)
  else: transtable += chr(min(max(n, 32),127))

class Data(StringLike):
    __constructor__ = "constructor_fromany"
    """Data objects behave as strings or buffers. It should be used to store
    binary data or file contents. It cannot be printed directly, but its
     contents can be accessed with the data() or textdata() methods

    It has a convert method, but the conversion engine is never used,
     convert() on a Data object is equivalent to cast()
    This behavior is NOT inherited if Data is Spyder-subclassed,
     in that case, convert() works normally

    It can be constructed from String or Data, but there is no parsing constructor
    In File objects, Data is considered to be a general format,
     convertible into and from everything
    """
    @staticmethod
    def typename():
      return "Data"
    def __eq__(self, other):
      if other is None: return False
      return self.data() == other.data()
    def __ne__(self,other): return not self.__eq__(other)  
    
    def __init__(self, s):
      return Data.__dict__[Data.__constructor__](self, s)
    def constructor_fromany(self, s):
      try:
        buf = s._buffer.raw
        self._buffer = s._buffer
      except AttributeError:
        try:
          buf = s.raw
          self._buffer = s
        except AttributeError:
          if len(s) and s[0] == s[-1]:
            s = String(s)
          try:
            s2 = bytes(s)     
          except TypeError:
            s2 = bytes(s, "utf-8")
          self._buffer = ctypes.create_string_buffer(s2)    
    def cast(self, target):
      return target(self._buffer.raw[:-1])
    def convert(self, target):
      """Data class conversion method
          Equivalent to a cast, unless Data is Spyder-subclassed"""
      if self.typename() == "Data":
        return self.cast(target)
      else: return spyder.core.convert(type(self),target, self, globals())
    def __str__(self):
      return String(repr(self))
    def __print__(self, spaces, mode):
      if mode == "str":
        return "<Data object of length %d>" % self.length()
      else:
        return repr(self)
    def __repr__(self):
      return repr(self._buffer.raw[:-1])
    def length(self):
      return len(self)
    def __len__(self):
      return max(len(self._buffer.raw) - 1, 0)
    def validate(self):
      """<empty>"""
      pass
    def data(self):
      """Returns the underlying data as binary C string"""
      return self._buffer.raw[:-1]
    def textdata(self):
      """Returns the underlying data as Python ASCII string:
          non-text characters are converted to whitespace"""
      s = str(self._buffer.raw[:-1])
      return s.translate(transtable)
    @classmethod
    def fromfile(c, filename, fastparse=False):
      return spyder.__types__["File"](filename, c).data()
    def tofile(self, filename):
      f = spyder.__types__["File"](filename, type(self), "w", self)
      f.save()
      f.close()
    def totempfile(self):
      from spyder.loader import file_exists
      import random, os,tempfile
      tempdir = tempfile.gettempdir()   
      while 1:
        tempnam = tempdir + os.sep + str(random.randrange(1,1000000000))+".tmp"
        if not file_exists(tempnam): break
      f = spyder.__types__["File"](tempnam, type(self), "w", self)
      f.temporary(True)
      f.save()
      return f
    def dict(self): 
      """Called by the dict function of Spyder classes
       that have Data members
      For internal use only"""
      return self.data()
      

spyder.__types__["Data"] = Data
spyder.__types__["ResourceData"] = _Resource(Data)
spyder.core.arrayfunc("Data", globals())
