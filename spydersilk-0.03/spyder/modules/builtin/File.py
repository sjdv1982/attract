# Copyright 2007, 2013 Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 

class File(Object):
    """File class: can read/save files and read URLs
    Permissions can be governed through various switches in system.fly
    Upon construction, the file <name> is read and its contents stored in <data>
    All File objects with the same <name> point to the same data
    """
    __constructor__ = "constructor_fromany"
    @staticmethod
    def typename(): return "File"
    def temporary(self, flag):
      """Set a "temporary" flag as signal for applications
      indicating that <data> should reside in memory,
      not on disk, and that <name> can be safely deleted
      """
      assert flag == True or flag == False
      self.__temporary = flag
    def is_temporary(self):
      """Query if the "temporary" flag was set"""
      return self.__temporary
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
      ret = {}
      ret["name"] = self.name
      ret["mode"] = self.mode
      ret["fileformat"] = self.fileformat().typename()
      ret["format"] = self.format().typename()
      return ret
    def __init__(self, *args, **kargs):
      return File.__dict__[File.__constructor__](self, *args, **kargs)
    def constructor_fromany(self, *args, **kargs):
      """Public general constructor,
          implementing the following specialized constructors:
        Keyword/value constructor:
         File(name, fileformat, mode="Read", data=None, format=None)
          <name>: file name or URL name
          <fileformat>: The fileformat of the data from <name>
            <fileformat> must be a Spyder type,
             not a Python type or a string constant
            if fileformat is Data, it is changed to the new type
             of <data> as soon as <data> is converted          
          <mode>: "Read"/"read"/"r" or "Write"/"write"/"w"
            URLS can never be opened in write mode
            Allowed values for local files and URLs can be controlled
             through switches in system.fly
          <data>: only allowed in write mode
            The contents of <name> are ignored, and the
            underlying object <data> is copied from this argument instead
          <format>: only allowed in write mode
            Sets the current format of <data>
        Copy constructor: File(<File>)
        Dict constructor: File({'name':..., ...})
        List constructor: File([name, ...])
        A parsing constructor is yet to be implemented, but not deemed necessary
         
      """ 
      try:
        if len(kargs) > 0 or len(args) != 1: raise TypeError
        self.__spydercopy__(args[0])
      except TypeError:
        try:
          if len(kargs) > 0 or len(args) != 1: raise TypeError
          self.__construct__(**args[0]) #dict
        except TypeError:                
          if len(kargs) == 0 and len(args) == 1: raise
          try:
            self.__construct__(*args, **kargs)
          except TypeError:
            if len(args) > 1 or len(kargs) > 0: raise
            self.__unpack__(args[0])
      self.__validate__()
    def __spydercopy__(self, file0):
      """Private copy constructor, for internal use only"""
      try:
        name = file0.name
        fileformat = file0.fileformat()
        mode = file0.mode
        temporary = file0.is_temporary()        
      except AttributeError:
        raise TypeError
      if name not in _resources: 
        file0.close()
        file0 = type(file0)(name, fileformat, mode)
      return self.__construct__(name = file0.name, mode=file0.mode, fileformat = file0.fileformat(), temporary=temporary)         
    def __construct__(self, *args, **kargs):
      """Private general constructor, for internal use only"""            
      name = None
      fileformat = None
      mode = None
      format = None
      data = None
      temporary = None
      if len(args) > 0:
        name = String(args[0])
      if len(args) > 1:
        fileformat = args[1]      
      if not fileformat is None:         
        if isinstance(fileformat, str): fileformat = spyder.__types__[String(fileformat)]
        assert type(fileformat) == type(int) and hasattr(fileformat, "typename")
      #it must be a Spyder type
      if len(args) > 2:
        mode0 = String(args[2])
        readcodes = ("r", "Read", "read")
        writecodes = ("w", "Write", "write")
        if mode0 not in readcodes and mode0 not in writecodes: raise TypeError("Unknown mode")
        if mode0 in writecodes: mode = "Write"
        else: mode = "Read"
      if len(args) > 3:
        data = args[3]  
      if len(args) > 4:
        format = args[4]
        if format != None:                 
          if isinstance(format, str): format = spyder.__types__[String(format)]
          assert type(format) == type(int) and hasattr(format, "typename")
      if len(args) > 4:
        temporary = args[5]
        assert temporary in (True, False)
            
      if "name" in kargs:
        if name != None: raise TypeError
        name = String(kargs["name"]) 
      if "fileformat" in kargs:
        if fileformat != None: raise TypeError
        fileformat = kargs["fileformat"]
        if fileformat != None:    
          if isinstance(fileformat, str): fileformat = spyder.__types__[String(fileformat)]          
          assert type(fileformat) == type(int) and hasattr(fileformat, "typename")
      if "mode" in kargs:        
        if mode != None: raise TypeError
        mode0 = String(kargs["mode"])
        readcodes = ("r", "Read", "read")
        writecodes = ("w", "Write", "write")
        if mode0 not in readcodes and mode0 not in writecodes: raise TypeError
        if mode0 in writecodes: mode = "Write"  
        else: mode = "Read"
      if "format" in kargs:
        if format is not None: raise TypeError
        format = kargs["format"]
        if format is not None:         
          if isinstance(format, str): format = spyder.__types__[String(format)]
          assert type(format) == type(int) and hasattr(format, "typename")         
      if "data" in kargs:
        if data != None: raise TypeError
        data = kargs["data"]
      if "temporary" in kargs:
        if temporary != None: raise TypeError
        temporary = kargs["temporary"]
        assert temporary in (True, False)                
      if mode == None: mode = "Read" #default is "Read"
      if temporary == None: temporary = False
      self.__temporary = temporary
      
      #print "name", name, "fileformat", fileformat, "mode", mode, "data", data, "format", format, "nodata", data == None

      if name == None: 
        raise TypeError("File must have a name")

      if mode == "Read" and data != None: 
        raise TypeError("Data must be undefined in Read mode")  
              
      name_loaded = False
      if name in _resources:
        name_loaded = True
         
      data_is_buffer = False
      data_loaded = False
      if data is not None:
        if format is not None and getattr(data,"typename", lambda: None)() != format.typename():
          data = format(data)
        if hasattr(data,"data") and hasattr(data,"typename") and not hasattr(data, "name"):
          data_is_buffer = True
        else:
          if hasattr(data,"name") and hasattr(data,"typename"):
            if data in _resources:
              data = data.name
              data_loaded = True
          else: data_is_buffer = True
            
      #print "name-loaded", name_loaded, "data-is-buffer", data_is_buffer, "data_loaded", data_loaded

      if name_loaded == True: #we are working with an existing File object
        self.name = name
        f = _resources[name]
        if data != None: #and mode is write, has been checked!
          if data_loaded == False and data_is_buffer == False and format == None: 
            raise TypeError
          loadfile(self.name, True, True) #check if we have write access, but do not re-read in
          if data_loaded == True: #copy from named File object <data>
            d = _resources[data]
            f = type(d)(d)
            fformat = type(d)
          elif data_is_buffer == True: #copy from buffer <data>, with specified format <format>
            fformat = type(data)
            f = fformat(data)       
          self.mode = mode
          self._fileformat = fileformat
          if fileformat == None:
            self._fileformat = fformat
        else: #no data, so basically nothing happens perhaps except a read=>write mode change
          if fileformat == None: raise TypeError
          self._fileformat = fileformat             
          self.mode = mode  
          if self.mode == "Write": 
            loadfile(self.name, True, True)
          f = _resources[name]
        #copy data and format
        self._refetodata = f
        _resources[name] = f
        if format != None: #convert to format <format>
          self.convert(format)
      else:   #we are creating a new File object
        self.name = name
        self.mode = mode
        if data is not None: #and mode is write, has been checked!
          loadfile(self.name, True, True) #only check if we have write access
          if data_loaded == False and format == None:
            if data_is_buffer == False:
              if fileformat == None: raise TypeError("You must supply a data format or a file format")
              else: format = fileformat                                                     
            else: 
              format = type(data)
              if isinstance(data, str): format = Data
          if data_loaded == True: #copy from named File object <data>
            d = _resources[data]
            if fileformat == None:  fileformat = type(d)
            self._fileformat = fileformat
            self._refetodata = d
            _resources[name] = d
            if format != None: self.convert(format)
          elif data_is_buffer == True: #copy from buffer <data>, with specified format <format>
            if isinstance(data, str):
              self._refetodata = Data(data)            
            else:
              self._refetodata = type(data)(data)            
            _resources[name] = self._refetodata 
            if fileformat == None: fileformat = format
            self._fileformat = fileformat
            if format != None: self.convert(format)
            else: self.convert(fileformat)
        else: #a new file without supplied data, just read the file in the specified format...
          self._fileformat = fileformat  
          if fileformat == None: raise TypeError("You must supply a file format")
          if format == None: format = fileformat
          if mode == "Read":
            d = loadfile(name, False, False)
            self._refetodata = fileformat(d) #to prevent the weakref dict to remove the data as long as this file exists
            _resources[name] = self._refetodata #just read it in
            
          else:
            loadfile(name, True, True) #check if we have write permission
            self._refetodata = fileformat("")
            _resources[name] = self._refetodata #initialize as empty buffer
            try: #try to read from file, but no problem if file does not (yet) exist
              d =  loadfile(name, False, False)
              self._refetodata = fileformat(d)
              _resources[name] = self._refetodata #just read it in
            except IOError:
              pass
          self.convert(format)          
    def __unpack__(self, s):
      """Private value/keyword/list constructor, for internal use only"""
      args,kargs = spyder.core.parse(s, "dict")
      if len(args) > 1:
        args[1]= spyder.__types__[args[1]]
      if len(args) > 3:
        args[3]= spyder.__types__[args[3]]
      if "fileformat" in kargs:
        kargs["fileformat"] = spyder.__types__[kargs["fileformat"]]
      if "format" in kargs:
        kargs["format"] = spyder.__types__[kargs["format"]]
      return self.__construct__(*args, **kargs)
    def data(self):
      """Returns <data>"""
      return _resources[self.name]
    def fileformat(self):
      """Returns <fileformat>, the original format of <data>"""
      return self._fileformat
    def format(self):
      """Returns <format>, the current format of <data>"""
      return type(_resources[self.name])
    def convert(self, newformat): 
      """Operates on the underlying <data>, not on the current object
      The underlying <data> is converted to "newformat"
      """
      if _resources[self.name].typename() != newformat.typename():
        if self._fileformat == Data:
          self._fileformat = newformat
        d = _resources[self.name].convert(newformat)
        self._refetodata = d
        _resources[self.name] = d      
      return _resources[self.name]
    def cast(self, newformat):
      """Operates on the underlying <data>, not on the current object
      The underlying <data> is cast to "newformat" """
      if _resources[self.name].typename() != newformat.typename():
        if self._fileformat == Data:
          self._fileformat = newformat
        d = _resources[self.name].cast(newformat)
        self._refetodata = d
        _resources[self.name] = d
      return _resources[self.name]     
    def length(self):
      """Returns the length of <data>"""
      return _resources[self.name].length()
    def __print__(self, spaces, mode):
      """Print handle for Spyder classes that embed File, for internal use only"""
      if mode == "str":
        return self.__str__(spaces)
      elif mode == "repr": 
        return self.__repr__(spaces)
    def __str__(self,spaces=0):
      """Triggers on str(self) and on print self
      Returns a string representation
      This representation can be parsed with spyder.core.parse
      """      
      sp = (spaces+2) * " "
      sp2 = spaces * " "
      s = "%s (\n%sname = \'%s\',\n%sfileformat = %s,\n%smode = \'%s\',\n%sformat = %s,\n%s)" %  (self.typename(), sp, self.name, sp, self.fileformat().typename(), sp, self.mode, sp, self.format().typename(), sp2)
      return s
    def __repr__(self,spaces=0):
      """Triggers on repr(self)
      Returns a parsable string representation and saves the current object"""
      s = self.__str__(spaces)
      self.save()
      return s
    def reload(self):
      d = loadfile(self.name, False, False)
      self._refetodata = self._fileformat(d) 
      _resources[self.name] = self._refetodata 
      self.convert(self.format())      
    def save(self):
      """
      If <mode> is "Write", saves <data> to disk
       but first, it converts <data> to <fileformat>
       unless <fileformat> is "Data"
      """       
      if self.mode == "Write":
        self.convert(self._fileformat)
        buf = self.data()
        savefile(self.name, buf)
    def close(self):
      """Unloads <data> from memory
      The current object and all other objects
       with the same <name> become invalid"""
      _resources.pop(self.name)
      self._refetodata = None      
    def delete(self):
      """Erases <name> from disk"""
      os.remove(self.name)
    @classmethod
    def fromfile(c, filename, fastparse=False):
      if fastparse: return c.fromdict(spyder.core.fastparse(filename, c)[1])
      else: return spyder.__types__["File"](filename, c).data()
    def tofile(self, filename):
      f = spyder.__types__["File"](filename, type(self), "w", self)
      f.save()
      f.close()
    def __eq__(self, other):
      try:
        return self.name == other.name and self.mode == other.mode
      except AttributeError:
        return False 
    def __ne__(self,other): return not self.__eq__(other)  
   
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
      
File.constructor_fromdict = File.constructor_fromany
spyder.__types__["File"] = File
spyder.core.arrayfunc("File",  globals())
