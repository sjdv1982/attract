import Spyder, spyder, sys
import traceback, inspect
import copy

from .resourcewrapper import embedded_resource_filename
from .model_status import model_status

_update_ok = "Model update OK"

def dict2list(d, full):
  if d is None: return None
  keys = set(d.keys())
  for k in keys: assert isinstance(k, int), k
  if full:
    ret = []    
    for n in range(max(keys)+1):
      if n in keys: ret.append(d[n])
      else: ret.append(None)
    return ret  
  l = len(keys)
  if len(keys) < max(keys) + 1:
    for n in range(len(keys)):
      if n not in keys: break         
    l = n
  if l == 0: 
    ret = None
  else: 
    ret = [d[n] for n in range(l)]    
  return ret

def determine_value(m, v):
  if m._is_resource:
    if isinstance(v, Spyder.File): v = (v.name, v.fileformat())
    elif isinstance(v, Spyder.Resource):
      if v.filename is not None:
        v = (v.filename, v.fileformat().typename())
    
    if isinstance(v, str) and v == embedded_resource_filename: return None, True
    if isinstance(v, tuple) and len(v) == 2 and \
      isinstance(v[0], str) and v[0] == embedded_resource_filename: 
      return None, True
  elif m._is_file:
    if isinstance(v, Spyder.File): v = (v.name, v.fileformat().typename())
    if v is None: return None, True
    if isinstance(v,str) and v=="": return None, True
    assert isinstance(v, tuple) and len(v) == 2 and isinstance(v[0], str), v
    filename = v[0]
    if isinstance(filename, str) and filename == embedded_resource_filename: return None, True
  return v, False
  
class model(object):
  _basic = False
  _is_file = False
  _is_resource = False
  _is_resource_file = False
  _is_valued = False
  _is_optional = False
  _requiredmembers = set()
  _async_children = set()
  _typedefaults = {}
  _constructionstate = "empty"
  _syncing = False
  _need_self_update = False
  
  def _start_syncing(self):
    assert not self._syncing
    if self._need_self_update: self._self_update()
    self._syncing = True

  def _stop_syncing(self):
    assert self._syncing
    if self._need_self_update: 
      self._self_update()
    self._syncing = False
  
  def __init__(self,type,parent=None,membername=None,value=None,is_default=False):    
    assert (parent is None) == (membername is None), (parent, membername)    
    self._toplevel = (parent is None)
    
    typetree = None
    ok = False
    try:
      if issubclass(type, Spyder.Object):
        self._type = type
        typename = type.typename()
        ok = True
    except TypeError:    
      pass   
    if not ok and isinstance(type, spyder.core.typetreeclass):       
      typetree = type
      typename = typetree.typename
      if hasattr(typetree, "type") and typetree.type is not None:        
        self._type = typetree.type
      elif typename is not None:
        self._type = getattr(Spyder, typename)
      else:
        self._type = None
      ok = True    
    if not ok and isinstance(type, str):  
      self._type = getattr(Spyder, type)
      typename = type
      ok = True
    if not ok:
      raise TypeError(type)
         
    if typename is not None and typename.startswith("Resource"):
      self._is_resource = True
      self._resourcetype = self._type
      typename = typename[len("Resource"):]
      self._type = getattr(Spyder, typename)
    
    if typetree is None and inspect.isclass(self._type):
      if hasattr(self._type, "_typetree"):
        typetree = self._type._typetree()              
      if issubclass(self._type, Spyder.File): self._is_file = True
        
    arraycount = 0
    if typename is not None:
      while typename.endswith("Array"):
        typename = typename[:-len("Array")]
        arraycount += 1
    self._typename = typename
    self._arraycount = arraycount
    self._parent = parent
    self._membername = membername
    
    if self._toplevel: self._path = ()
    else: self._path = self._parent._path + (self._membername,)
    
    self.__dif = None
    self.__ddif = None
    self.__exc = None   
        
    self.__passup = False
    self.__passdown = False
    self._value = None
    self._truevalue = None    
    
    self._children = {}
    self._async_children = set()
    if self._arraycount:      
      pass #should be taken care of by self._set, in just a moment
    elif typetree is not None and typetree.members is not None and len(typetree.members):
      self._typedefaults = {}
      self._requiredmembers = set()      
      for membername, mtypetree in typetree.members:
        mdefault = mtypetree.default        
        self._typedefaults[membername] = mdefault
        mvalue = mdefault
        if value is not None:
          mvalue = getattr(value, membername) 
        mtypename = None
        if mtypetree.typename is None:
          child = model(mtypetree, self, membername, mvalue,mtypetree.is_default)
        else:  
          mtypename = mtypetree.typename + mtypetree.arraycount * "Array"  
          if mtypetree.is_resource: mtypename = "Resource" + mtypename
          if hasattr(mtypetree,"typemapped") and mtypetree.typemapped: mtypename = mtypetree 
          child = model(mtypename, self, membername, mvalue,mtypetree.is_default)
        if not child._is_optional: self._requiredmembers.add(membername)
        self._children[membername] = child        
    else: #basic model
      self._basic = True
      
    #now, determine our optionality (_is_optional)
    self._is_default = is_default
    if is_default:
      self._is_optional = True
    elif self._arraycount:
      self._is_optional = True
    elif self._complete():
      if self._type is not None:
        v = None
        self._constructionstate = "full"
        try:
          v = self._type.fromdict(self._std())        
        except:  
          self._constructionstate = "failed"
          self.__exc = sys.exc_info()
        if self._toplevel:
          self._truevalue = v
        else:
          self._is_optional = True  
      else:
        self._constructionstate = "full"
    
    self._listeners = []
    self._set(value, _toplevelcall=False)
    self.__passdown = True    
    if self._constructionstate == "empty":
      for childname in self._children:
        child = self._children[childname]
        if child._constructionstate == "full":
          self._update_child(childname)
    self.__passup = True
    
  def _listen(self, listener):
    self._listeners.append(listener)
  def _unlisten(self, listener):
    if listener in self._listeners: 
      self._listeners.remove(listener)    
      return True
    else: return False
    
  def __update_dif(self, childname, childvalue, childdif, childarraycount):
    if not self._arraycount:
      if childvalue == self._typedefaults[childname]: childvalue = None
    if childvalue is None and self.__dif is None: return False
    
    change = True
    dif = self.__dif
    ddif = self.__ddif
    if dif is None:
      dif = {}
    if ddif is None:  
      ddif = {}
    if childvalue is None:      
      if childname not in dif: 
        change = False
      else:
        dif.pop(childname)
        self._child_sync(childname)
        if not len(dif): 
          dif = None        
      if len(ddif):    
        if childname in ddif: change = True
        ddif.pop(childname, None)
        if not len(ddif):   
          ddif = None     
          if len(self._requiredmembers): 
            self._constructionstate = "empty"
    else:  
      updif = False
      if childname in dif and dif[childname] == childvalue: 
        change = False
      else:  
        updif = True
        dif[childname] = type(childvalue)(childvalue)
      if childdif is not None:
        ddif[childname] = copy.copy(childdif)
        change = True
        self._constructionstate = "partial"
      elif updif:  
        ddif[childname] = {}
      elif childarraycount:
        ddif[childname] = type(childvalue)(childvalue)
      else:
        if not len(ddif): ddif = None
    self.__dif = dif
    self.__ddif = ddif
    return change
   
  def _update_children(self): 
    if not self.__passdown: return
    #self._async_children.clear()      
    keys = set(self._children.keys())
    if self._arraycount and self._truevalue is not None:
      keys = keys.union(set(range(len(self._truevalue))))
    for childname in keys:
      mvalue = self._get_child_value(childname)
      child = self._get_child(childname)
      child._disable_passup()
      child._set(mvalue, _toplevelcall=False)
      child._enable_passup()
    
  def _update_child(self, childname):        
    #This function is called by the child, to notify the parent        
    child = self._children[childname]        
    change = self.__update_dif(childname, child._value, child._dif(), child._arraycount)
    if change:
      if self._syncing: 
        self._need_self_update = True
      else:
        self._self_update()
  
  def _self_update(self):    
    self._need_self_update = False
    #build value from dif
    self.__exc = None
        
    if not self._complete(): 
      if self.__ddif is not None: return
      self._constructionstate = "empty"
      try:      
        if self._arraycount:
          self._value = self._type([])
        else:
          self._value = self._type(self._std())
      except:
        self._value = None
    else:
      self._constructionstate = "failed"     
      if self._type is None:
        d = copy.copy(self.__dif)
      else:
        d = copy.copy(self.__ddif)
      if self._arraycount:
        d = dict2list(d, full=False)      
      try:
        if self._arraycount or len(self._children):
          for childname, child in self._children.items():
            if child._arraycount and not child._is_default and childname not in d:
              d[childname] = []
          if d is None: 
            if not self._toplevel and (self._parent._arraycount or self._parent._typedefaults[self._membername] is None):
              value = None
            elif self._arraycount:
              value = self._type([])
            elif self._type is None:
              value = None
            else:
              value = self._type(d)
          else: 
            for childname, child in self._children.items():
              if child._arraycount and childname in self._requiredmembers and childname not in d: d[childname] = []
            if self._type is None:
              value = d
            else:
              value = self._type.fromdict(d)
        else:
          value = self._type(d)
        self._truevalue = value
        if not self._toplevel and not self._parent._arraycount:
          if self._parent._typedefaults[self._membername] == value:
            value = None
        self._value = value      
        self._constructionstate = "full" 
        if not len(self._requiredmembers) and self.__ddif is None: 
          self._constructionstate = "empty"
      except:      
        self.__exc = sys.exc_info()
        if not self._toplevel:
          self._parent._child_async(self._membername)
        return

    if self._type is not None:
      self._update_children()
    
    #Are we valued or unvalued?
    if self.__ddif is not None or self._toplevel: 
      self._is_valued = True    
    else:
      self._is_valued = False
    
    if not self._toplevel and not len(self._async_children):
      self._parent._child_sync(self._membername)      
    for listener in self._listeners:      
      listener(self._truevalue)          
    if self.__passup and not self._toplevel:
      self._parent._update_child(self._membername)
     
  def _child_async(self, childname):
    child = self._children[childname]
    #if not child._is_valued: return
    if not len(self._async_children) and not self._toplevel:
      self._parent._child_async(self._membername)
    self._async_children.add(childname)
    
  def _child_sync(self, childname):
    if childname in self._async_children:
      self._async_children.remove(childname)      
    if not len(self._async_children) and not self._toplevel:
      self._parent._child_sync(self._membername)
      
  def _get_child(self, child):
    if isinstance(child, int):
      if child not in self._children:
        for n in range(child+1):
          if n not in self._children: 
            subtypename = self._typename + (self._arraycount - 1) * "Array"
            self._children[n] = model(subtypename, self, n, None, False)
    if child not in self._children: raise AttributeError(child)
    return self._children[child]
  
  def _get_child_value(self, child):
    if self._truevalue is None: return None
    if self._arraycount: 
      ret = None
      if child < len(self._truevalue):
        ret = self._truevalue[child]
    elif self._type is None:
      if child not in self._truevalue: return None
      ret = self._truevalue[child]
    else: 
      ret = getattr(self._truevalue, child)
    return ret
  
  def append(self, item):
    assert self._arraycount > 0
    index = len(self._children)
    return self._get_child(index)._set(item)
  
  def _insert(self, index):
    assert self._arraycount > 0
    if index > len(self._children): return
    for n in reversed(range(index, len(self._children))):
      if n in self._children:
        child = self._children.pop(n)
        child._membername = n+1
        self._children[n+1] = child
      if self.__dif is not None and n in self.__dif:
        child = self.__dif.pop(n)
        self.__dif[n+1] = child
      if self.__ddif is not None and n in self.__ddif:
        child = self.__ddif.pop(n)
        self.__ddif[n+1] = child
    subtypename = self._typename + (self._arraycount - 1) * "Array"
    self._children[index] = model(subtypename, self, index, None, False)

  def _delete(self, index):
    assert self._arraycount > 0
    if index >= len(self._children): return False    
    l = len(self._children) - 1
    maxchild = 0
    for n in range(l+1):
      if n in self._children:
        child = self._children[n]
        if child._constructionstate != "empty": maxchild = n
    for n in range(index, l-1):
      self._children.pop(n, None)
      if self.__dif is not None:
        self.__dif.pop(n, None)
      if self.__ddif is not None:
        self.__ddif.pop(n, None)
      if (n+1) in self._children:
        child = self._children.pop(n+1)        
        child._membername = n
        self._children[n] = child
      if self.__dif is not None and (n+1) in self.__dif:
        child = self.__dif.pop(n+1)
        self.__dif[n] = child
      if self.__ddif is not None and (n+1) in self.__ddif:
        child = self.__ddif.pop(n+1)
        self.__ddif[n] = child
    if self.__dif is not None and not len(self.__dif): self.__dif = None
    if self.__ddif is not None and not len(self.__ddif): self.__ddif = None
    subtypename = self._typename + (self._arraycount - 1) * "Array"
    self._children[l] = model(subtypename, self, l, None, False)
    
    ls = len(self)    
    if self.__ddif is None or ls < maxchild: 
      return False
    
    value = self._type(dict2list(self.__dif,full=False))
    for n in range(index, l):
      child._disable_passup()
    self._set(value)  
    for n in range(index, l):      
      child._enable_passup()
    return True
    
  def __getitem__(self, index):
    assert isinstance(index, int)
    assert self._arraycount > 0
    assert index >= 0 #no negative indices
    return self._get_child(index)  
  def __setitem__(self, index, value):
    assert isinstance(index, int)
    assert self._arraycount > 0
    assert index >= 0 #no negative indices
    return self._get_child(index)._set(value)
    if value is None:
      for n in range(index,-1,-1):
        if n not in self._children: continue
        if self._children[n] is not None: continue
        self._children.pop(n)
        if self._dif is not None: self._dif = self._dif[:index]
        if self._ddif is not None: self._ddif = self._ddif[:index]
  def __contains__(self, key):
    return key in self._children
  def __delitem__(self, key):
    if key in self._children: 
      self._children._set(None)
      del self._children[key]
  def __iter__(self):
    assert self._arraycount > 0
    for n in range(max(self._children.keys())+1):   
      child = self._children[n]        
      if child._constructionstate == "empty": break 
      yield child
  def __len__(self):
    assert self._arraycount > 0
    l = 0
    for n in range(max(self._children.keys())+1):   
      child = self._children[n]        
      if child._constructionstate == "empty": break 
      l += 1
    return l  
  def __getattr__(self, attr):  
    if attr.startswith("_"): raise AttributeError(attr)
    return self._get_child(attr)
  def __setattr__(self, attr, value):
    if attr.startswith("_"): 
      self.__dict__[attr] = value
    else:  
      self._get_child(attr)._set(value, _toplevelcall = False)      
      return self._status()

  def _set(self, value, _toplevelcall=True): 
    need_update = False
    if isinstance(value, model): value = value._get()
    #print("model._set", self._path, str(value)[:30])    
    value, cancel = determine_value(self, value)
    if cancel: return self._status()
    if value == self._value: return _update_ok
    
    is_resource_file = False
    if self._is_resource:
      if isinstance(value, tuple) and len(value) == 2:
        is_resource_file = True
        if value[0] == None: value = None
        elif value[0] == "": value = None        
      
    self.__exc = None
    try:      
      if self._basic:
        
        #1: basic model, can be file or resource        
        if value is not None: self._constructionstate = "failed"
        if self._is_file or is_resource_file:
          filename, fileclassname = value                  
        if self._is_file:  
          fileclass = getattr(Spyder, fileclassname)        
          v = Spyder.File(filename,fileclass,"r")
        elif self._is_resource:
          if is_resource_file:
            v = self._resourcetype(filename=filename)
            self._is_resource_file = True
          else:
            v = self._resourcetype(value)
        elif inspect.isclass(self._type) and issubclass(self._type, Spyder.Bool):
          v = True if value else False
        else:
          v = self._type(value)      
        value = v  
        self._constructionstate = "full" if value is not None else "empty"
      elif self._is_resource:
        #2: non-basic resource model
        self._is_resource_file = False      
        if is_resource_file:
          filename,fileclassname = value
          fileclass = getattr(Spyder, fileclassname)
          fileobj = Spyder.File(filename,fileclass,"r")
          v0 = self._resourcetype(fileobj)
          self._is_resource_file = True
          self._constructionstate = "full"
        else:
          v0 = self._resourcetype(value)    
          if v0.filename is not None:
            self._is_resource_file = True
          self._constructionstate = "full" if value is not None and value != {} else "empty"  
        if self._is_resource_file: 
          value = v0.file()
          self._constructionstate = "full"
        else:
          v = v0.data()
          self._constructionstate = "full" 
          if not self._toplevel and value is None: self._constructionstate = "empty"
        value = v  
      elif value is None and self._arraycount:
        value = self._type([])
        self._constructionstate = "empty"
      elif self._type is None:
        if value is not None:
          assert isinstance(value, dict)
          for childname in self._children:
            if childname not in value: 
              self._constructionstate = "partial"
              return
          self._value = {}
          for childname in self._children:
            child = self._children[childname]
            child._disable_passup()
            child._set(value[childname], _toplevelcall=False)
            child._enable_passup()          
            self._value[childname] = child._value
          self._truevalue = self._value  
          self.__dif = copy.copy(self._value)  
          self.__ddif = copy.copy(self._value)  
          self._constructionstate = "full"  
          self._is_valued = True
          return _update_ok    
        else:
          self.__dif = None
          self.__ddif = None
          self._constructionstate = "empty"
          if self._complete(): self._constructionstate = "full"
          self._is_valued = False          
      else:    
        #non-basic non-resource model
        value0 = value
        value = self._type(value)
        self._constructionstate = "full"
        if not self._toplevel and (value0 is None or value0 == {}): self._constructionstate = "empty"
    except:
      self.__exc = sys.exc_info()
      #return, unless we just cleared the value
      if value is None and self._value is not None:
        self._constructionstate = "empty"        
      else:
        if not self._toplevel:
          self._parent._child_async(self._membername)
        if _toplevelcall: return self._status()
        return
    
    if not self._toplevel and not len(self._async_children): 
      self._parent._child_sync(self._membername)
    if self._truevalue != value: 
      need_update = True
    self._truevalue = value  
    if not self._toplevel and not self._parent._arraycount:
      if self._parent._typedefaults[self._membername] == value:
        value = None
    self._value = value
    self._update_children()
    
    self.__dif = None
    self.__ddif = None
    c = self._constructionstate
    if self._truevalue is not None:
      if self._arraycount:
        for n in range(len(self._truevalue)):          
          self.__update_dif(n, self._truevalue[n], self._truevalue[n].dict(), self._arraycount-1)
        self._constructionstate = c
      elif len(self._children) > 0 and not self._is_resource_file:
        for childname in self._children:
          mvalue = self._get_child_value(childname)
          child = self._get_child(childname)
          self.__update_dif(childname, mvalue, child._dif(), child._arraycount)
        self._constructionstate = c  
      elif self._is_file or self._is_resource:
        self.__dif = self._value.dict()
        self.__ddif = self._value.dict()        
      else:
        self.__dif = self._value
        self.__ddif = self._value
        
    #Are we valued or unvalued?
    if self.__exc is None:
      if self.__ddif is not None or self._toplevel: 
        self._is_valued = True    
      else:
        self._is_valued = False        
    else:
      self._is_valued = False
      
    for listener in self._listeners:
      listener(self._truevalue)        
    if self.__passup and not self._toplevel and need_update:
      self._parent._update_child(self._membername)
      
    return _update_ok
    
  def _get(self): 
    if self._type is None:
      if self._truevalue is None: return None
      ret = {}
      for k,v in self._truevalue.items():        
        child = self._children[k]
        vv = v
        if not child._basic and child._type is not None:
          vv = v.dict()
        ret[k] = vv
      return ret
    else:  
      if not self._is_valued: return None
      return self._truevalue
      
  def _dif(self, full=False): 
    if self._arraycount:
      return dict2list(self.__ddif, full=full)
    else:      
      ret = self.__ddif      
      if full:
        if not isinstance(ret, dict): return ret
        ret2 = {}
        for k in ret:
          ret2[k] = self._get_child(k)._dif(full=True)
        if not len(ret2): return None
        return ret2
      else:
        return ret
  def __str__(self):
    return str(self._get())
  def __repr__(self):
    return repr(self._get())
    
  def _enable_passup(self): self.__passup = True
  def _disable_passup(self): self.__passup = False
  def _enable_passdown(self): self.__passdown = True
  def _disable_passdown(self): self.__passdown = False
  
  def _std(self):
    std = {}
    for childname, child in self._children.items():
      if child._arraycount: std[childname] = []
    return std  
  
  def _complete(self):
    if self._basic: return False
    complete = True
    if len(self._requiredmembers):
      complete = False
      defined = set()
      if self.__dif is not None: defined = self.__dif.keys()
      if not self._requiredmembers.difference(defined): complete = True #all members are defined
    return complete
    
  def _status(self, controller=None):
    return model_status(self, controller)
    
  def _clear(self):
    self.__dif = None
    self.__ddif = None
    self.__exc = None           
    self._value = None
    v = None
    try:
      v = self._type(v)
    except:
      pass
    self._truevalue = v
    self._is_resource_file = False
    self._is_valued = False
    self._async_children = set()
    for child in self._children.values(): child._clear()
  
  def fromfile(self, filename, fastparse=False):
    assert self._type is not None
    obj = self._type.fromfile(filename, fastparse=fastparse)
    self._set(obj)
    
  def tofile(self, filename):
    assert self._type is not None
    self._get().tofile(filename)    