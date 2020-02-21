# Copyright 2013, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

import os, weakref

def loadfile(filename, checkwrite, checkonly):
  if checkwrite == True:
    assert checkonly == True
  filename2 = os.path.split(filename)
  if filename2[0] == "":
    if checkwrite:
      func = writeCurrDir
    else:
      func = readCurrDir
  elif filename.find("://") > -1:
    if checkwrite:
      func = writeRemote
    else:
      func = readRemote
  else:
    if checkwrite:
      func = writeLocal
    else:
      func = readLocal
  if checkonly == True:
    if func(filename, True) != True:
      if checkwrite == True:
        raise Exception("Unable to write to file %s" % filename)
      else:
        raise Exception("Unable to read from file %s" % filename)
    return 
  else:
    f = func(filename, False)
    if f == None:
      raise Exception("Unable to read from file %s" % filename)
    try:
      ret = f.read()
      #print "Successfully opened file %s" % (filename)
      return ret
    finally:
      f.close()

def savefile(filename, buf):
  filename2 = os.path.split(filename)
  if filename2[0] == "":
    func = writeCurrDir
  elif filename.find("://") > -1:
    func = writeRemote
  else:
    func = writeLocal
  f = func(filename, False)
  if f == None:
    raise Exception("Unable to write to file %s" % filename)
  try:
    buf._buffer
    s = buf.data()
  except:
    try:
      buf.typename()
      s = repr(buf).encode('UTF-8')
    except:
      s = str(buf).encode('UTF-8')
  try:
    f.write(s)
    #print "Successfully saved file %s" % filename
  finally:
    f.flush()
    f.close()
  
_resources = weakref.WeakValueDictionary()

class Resource(Object):
  __constructor__ = "constructor_fromany"
  @classmethod
  def typename(cls): return "Resource" + cls.__resourceformat__.typename()
  
  @property
  def name(self): return self.filename
  @name.setter
  def name(self, name): self.filename = name 
  
  def __validate__(self):
    """<empty>"""
    pass
  def validate(self):
    """<empty>"""
    pass  
  def length(self):
    """Redirects to __len__"""
    return len(self)       
  def dict(self):
    """Returns a Python dict representation of the current object"""
    if self._data is not None and self.filename is not None:
      raise ValueError("Cannot represent %s: can be interpreted as filename or embedded object"  % self.typename())
    if self._data is not None:
      self._sync()      
    ret = dict()      
    if self.filename is not None: #file-like representation
      ret["name"] = self.filename
      ret["fileformat"] = self.__resourceformat__.typename()
      ret["mode"] = "Read"
      ret["format"] = self.__resourceformat__.typename()
    else: #embedded representation
      return self._data.dict()
    return ret  
  def __init__(self, *args, **kargs):
    cls = self.__class__
    return getattr(cls, cls.__constructor__)(self, *args, **kargs)
      
  def __construct__(self, data, filename, fileformat):
    """General private constructor, for internal use only"""
    if filename is None and data is None: 
      data = self.__resourceformat__()
    if filename is None and fileformat is not None: raise TypeError                  
    if data is not None:
      if isinstance(data, self.__resourceformat__):
        assert filename is None or filename not in _resources, filename
        data = self.__resourceformat__(data)        
      elif isinstance(data, File):
        if data.mode == "Write":
          data = data.data()
        else:
          filename = data.name
          fileformat = data.fileformat().typename()
          data = None
      elif isinstance(data, Resource):
        filename = data.filename
        fileformat = data.fileformat().typename()
        if filename is None: data = data._data 
        else: data = None        
      else:
        assert filename is None or filename not in _resources, filename
        data = self.__resourceformat__(data)
    self.filename = None
    if filename is not None: self.filename = String(filename)
    if fileformat is not None:
      if fileformat not in (self.__resourceformat__, self.__resourceformat__.typename()): raise ValueError(fileformat, self.__resourceformat__)        
    self._data = None
    if data is not None:
      self._data = data
      if filename is not None: 
        _resources[filename] = data        
        self._dataref = data #just to keep it alive
    self.__validate__()
    
  def constructor_fromdict(self, dic):
    """Dict constructor"""
    if isinstance(dic, str) or (spyder.python3 and isinstance(dic, bytes)): 
      return self.__construct__(String(dic), None, None)
    filename = None
    data = None
    fileformat = None  
    if "filename" in dic: filename = dic["filename"]
    elif "name" in dic: 
      if "mode" in dic and "data" in dic and dic["mode"] in ("w", "Write"):
        return self.__construct__(dic["data"], None, None)
      filename = dic["name"]      
    elif "data" in dic:
      return self.__construct__(dic["data"], None, None)
    else:
      return self.__construct__(dic, None, None)
    if "fileformat" in dic: fileformat = dic["fileformat"]
    return self.__construct__(None, filename, fileformat)
    
  def constructor_fromany(self, *args, **kargs):
    if len(kargs) == 0 and len(args) == 1 and isinstance(args, dict):
      return self.constructor_fromdict(args[0])
    filename = None
    data = None
    fileformat = None      
    if len(args) > 1:
      raise TypeError("Resource can have only one unnamed argument")
    if "filename" in kargs: filename = kargs["filename"]
    elif "name" in kargs: filename = kargs["name"]
    if "data" in kargs: data = kargs["data"]
    if "fileformat" in kargs: fileformat = kargs["fileformat"]
    if args:  
      if isinstance(args[0], str) and args[0].find("\n") == -1 and args[0].find(" ") == -1:
        filename = args[0]
      else:
        data = args[0]
    self.__construct__(data, filename, fileformat)
    
  def data(self):
    """Returns <data>"""
    self._sync()
    return self._data
  def fileformat(self):
    return self.__resourceformat__
  def format(self):
    return self.__resourceformat__      
  def convert(self, newformat): 
    """Operates on the underlying <data>, not on the current object
    The underlying <data> is converted to "newformat" and a copy is returned
    """
    d = self.data()
    return d.convert(newformat)
  def cast(self, newformat):
    d = self.data()
    return d.cast(newformat)
  def length(self):
    """Returns the length of <data>"""
    return self.data().length()
  def __print__(self, spaces, mode):
    """Print handle for Spyder classes that embed Resource, for internal use only"""
    if mode == "str":
      return self.__str__(spaces)
    elif mode == "repr": 
      return self.__repr__(spaces)
  def __str__(self,spaces=0,mode="str"):
    """Triggers on str(self) and on print self
    Returns a string representation
    This representation can be parsed with spyder.core.parse
    """
    if self._data is not None and self.filename is not None:
      raise ValueError("Cannot represent %s: can be interpreted as filename or embedded object" % self.typename())
    if self._data is not None:
      self._sync()      
    if self.filename is not None: #file-like representation
      s0 = "File"
      s1 = self.filename
      s2 = self.__resourceformat__.typename()
      s3 = "Read"
      s4 = self.__resourceformat__.typename()
      sp = (spaces+2) * " "
      sp2 = spaces * " "
      s = "%s (\n%sname = \'%s\',\n%sfileformat = %s,\n%smode = \'%s\',\n%sformat = %s,\n%s)" %  (s0, sp, s1, sp, s2, sp, s3, sp, s4, sp2)
      return s
    else: #embedded representation
      return self._data.__print__(spaces, mode)
    
  def __repr__(self,spaces=0):
    """Triggers on repr(self)
    Returns a parsable string representation and saves the current object"""
    s = self.__str__(spaces,mode="repr")
    if self.filename is not None: self.save()
    return s

  def file(self):
    assert self.filename is not None
    f = self.__resourceformat__
    return Spyder.File(name=self.filename,fileformat=f,format=f,mode="Read")
  
  def dict(self):
    if self._data is not None: 
      return self._data.dict()
    elif self.filename is not None:
      return self.file().dict()
    else: 
      return None

  def copy(self):
    return type(self)(self)
      
  def _sync(self):
    """
    If _resources exists, synchronizes _resources onto _data
    else, synchronizes _data onto _resources
    if neither exists, load _resources
    """
    if self.filename is not None and self.filename in _resources:  
      d = _resources[self.filename]
      if d is not self._data:
        if d.typename() != self.__resourceformat__.typename():
          d = d.convert(self.__resourceformat__)
          _resources[self.filename] = d
        self._data = d
        self._dataref = d #just to keep it alive
    elif self._data is not None:
      if self.filename is not None: 
        _resources[self.filename] = self._data
        self._dataref = self._data #just to keep it alive
    else:
      if self.filename is None: raise ValueError
      self.load() 
        
  def set(self, data):
    data = self.__resourceformat__(data)
    self._data = data
    if self.filename is not None:
      d = self.__resourceformat__(data)
      _resources[self.filename] = d
      self._dataref = d #just to keep it alive
  def load(self):
    if self.filename is None: raise ValueError
    data = self.__resourceformat__.fromfile(self.filename)
    if self.filename not in _resources or _resources[self.filename] != data: 
      _resources[self.filename] = data
      self._dataref = data #just to keep it alive
    self._data = data
  def save(self):
    if self.filename is None: raise ValueError
    if self._data is not None:
      raise ValueError("Cannot represent %s: can be interpreted as filename or embedded object"  % self.typename())    
    if self.filename not in _resources: return
    data = _resources[self.filename]
    if data.typename() != self.__resourceformat__.typename():
      data = data.convert(self.__resourceformat__)      
      _resources[self.filename] = data
      self._dataref = data #just to keep it alive
    savefile(self.filename, data)
      
  def delete(self):
    """Erases <name> from disk"""
    os.remove(self.filename)
  @classmethod
  def fromfile(c, filename, fastparse=False):
    return spyder.__types__["File"](filename, c).data()
  def tofile(self, filename):    
    data = self.data()
    data.tofile(filename)      
  def __eq__(self, other):
    try:
      return self.format() == other.format() and self._data == other._data and self.filename == other.filename
    except (TypeError,AttributeError,ValueError):
      return False  
  def __ne__(self,other): return not self.__eq__(other)
  
  def embed(self):
    """
    Interprets the resource as an embedded data object
    Eliminates any filename representation
    Loads the file if necessary
    """
    if self.filename is not None:
      self._sync()
      self.filename = None
      self._dataref = None
  def link(self, filename = None):
    """
    Interprets the resource as a filename
    Eliminates any embedded representation
    Does not save or load the file!
    """      
    if self._data is None: 
      if filename is not None and filename != self.filename: 
        self.filename = filename
        self.load()
        self._data = None
      return
    
    if filename is not None: self.filename = filename
    else: assert self.filename is not None
    _resources[self.filename] = self._data
    self._dataref = self._data #just to keep it alive
    self._data = None
    

  @classmethod
  def _form(cls):
    """
     Returns the spyderform of the class
    """
    return spyder.core.spyderforms[cls.typename()]
  @classmethod
  def _typetree(cls): 
    """ 
      Returns the typetree of the Spyder class
    """          
    return spyder.core.get_typetreedict(cls.typename())
    
    
def _Resource(spyderclass):
  import functools
  c = spyderclass
  if isinstance(c, str): c = spyder.__types__[c]
  typename = "Resource"+c.typename()
  cls = type(typename, (Resource,), dict(__resourceformat__=c))
  cls.fromdict = functools.partial(spyder.__constructor, 
    "constructor_fromdict",
    cls,
    "constructor_fromany",
  )
  spyder.__types__[typename] =  cls
  spyder.core.spyderform(cls._typetree())   
  return cls
  
spyder.__types__["_Resource"] = _Resource
spyder.__types__["Resource"] = Resource
