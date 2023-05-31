import Spyder
import cgi as stdcgi

def subdic(dic, name):
  #Gets a CGI subdict
  if name in dic:
    return True, dic[name]
  ret = {}
  name2 = name + "-"
  for k in dic:
    if k.startswith(name2):
      ret[k[len(name2):]] = dic[k]
  if not len(ret): return None, None
  return False, ret

def tostr(v):
  #Converts a default value to a string representation
  # or: converts a CGI empty string representation to None
  if isinstance(v, str): 
    if v == "": return None
    if v == "True": return True
    if v == "False": return False    
    return v
  if v is None: return None    
  if isinstance(v, bool):
    if v == True: return True
    if v == False: return False
  
  if isinstance(v, Spyder.File): return v
  try:
    if len(v) == 0: return None
  except (TypeError, AttributeError):
    pass   
  return tostr(str(v))

def list_tostr(lis):
  #Converts a default value list (obtained with SpyderType.list()/.dict())
  # to a string representation list
  ret = []
  for m in lis:
    if isinstance(m, dict):
      mm = dic_tostr(m)
    elif isinstance(m,list) or isinstance(m, tuple):
      mm = list_tostr(m)
    else:
      mm = tostr(m)
    ret.append(mm)
  if ret == [None] * len(lis): ret = None
  return ret
  
def dic_tostr(dic):
  #Converts a default value dict (obtained with SpyderType.dict())
  # to a string representation dict
  if isinstance(dic,list) or isinstance(dic, tuple):
    return list_tostr(dic)
  ret = {}
  for k in dic:
    m = dic[k]
    if isinstance(m, dict):
      mm = dic_tostr(m)
    elif isinstance(m,list) or isinstance(m, tuple):
      mm = list_tostr(m)
    else:
      mm = tostr(m)
    ret[k] = mm
  return ret
        

def _filter_delta(dic, spyderform,  typetree, default):
  """
  Filters a reduced CGI dict by setting all empty entries to None
  Returns:
    v: The reduced CGI dict
  """  
    
  if default is None: 
    try:
      default = spyderform.default
    except AttributeError:
      pass

  cls = getattr(Spyder, spyderform.typename)
   
  typ = None 
  if hasattr(spyderform, "type"):
    t = spyderform.type
    if t in ("none", None): return None
    typ = t
  if typ is None:
    if issubclass(cls, getattr(Spyder, "File")):
      typ = "file"
    elif issubclass(cls, getattr(Spyder, "Bool")):
      typ = "checkbox"
    elif issubclass(cls, getattr(Spyder, "Float")) or issubclass(cls, getattr(Spyder, "Integer")):  
      typ = "number"
    elif issubclass(cls, getattr(Spyder, "String")):  
      typ = "text"
  if typ is not None: 
    return dic 
   
  memnames = spyderform.get_membernames()
  if spyderform.arraycount > 0:
    if dic is None: return None
    ret = []
    length = max(spyderform.length, len(dic))
    for n in range(length):
      sub = None
      if dic is not None and len(dic) > n:
        sub = dic[n]
      mdefault = None
      if default is not None and len(default) > n: mdefault = default[n]        
      if sub is None:
        #absolutely no CGI entries for this index... must be non-existing in the form
        v = None
      else: 
        sf = spyderform[n] if n < spyderform.length else spyderform[0]
        v = _filter_delta(sub, sf, sf._typetree, mdefault)
      ret.append(v)
    pos = len(ret) - 1
    while pos >= 0 and ret[pos] == None: pos -= 1
    ret = ret[:pos+1]
    if ret == [None] * len(ret): ret = None
    v = ret
  elif not memnames:
    return dic
  else: #composite Spyder object 
    if dic is None: return None
    v = {}
    for membername in spyderform.get_membernames():
      if membername not in dic: continue
      member = getattr(spyderform, membername)
      sub = dic[membername]
      mtt = [e[1] for e in typetree.members if e[0] == membername]
      if len(mtt) == 0: continue #"virtual" form member
      if len(mtt) > 1: raise Exception(membername)
      mtt = mtt[0]
      mdefault = None
      if default is not None: mdefault = getattr(default, membername)
      vv = _filter_delta(sub, member, mtt, mdefault)
      if vv is not None: v[membername] = vv
    if not len(v): v = None
  
  return v
        
def _reduce_cgi(dic, spyderform, prename, is_value, typetree, default):
  """
  Reduces (a part of) a CGI dict
  Returns:
    v: The reduced CGI dict
  """  
  if default is None: 
    try:
      default = spyderform.default
    except AttributeError:
      pass

  prenam = ""
  if prename is not None: prenam = prename + "-"
          
  cls = getattr(Spyder, spyderform.typename)
    
  typ = None
  if hasattr(spyderform, "type"):
    t = spyderform.type
    if t in ("none", None): return None
    typ = t
  if typ is None:
    if issubclass(cls, getattr(Spyder, "File")):
      typ = "file"
    elif issubclass(cls, getattr(Spyder, "Bool")):
      typ = "checkbox"
    elif issubclass(cls, getattr(Spyder, "Float")) or issubclass(cls, getattr(Spyder, "Integer")):  
      typ = "number"
    elif issubclass(cls, getattr(Spyder, "String")):  
      typ = "text"
      
  memnames = spyderform.get_membernames()
  if typ == "file" and spyderform.is_resource:
    formdefault = None #a filled-in file is always a form modification
    if dic is None or dic == "":
      v = None
    else:  
      v = {"data": dic}
  elif memnames is None or not len(memnames) or typ is not None:
    if spyderform.arraycount > 0:
      assert is_value is not True, (prename, spyderform.typename)
      assert hasattr(spyderform, "length"), (prename, spyderform.typename)
      if dic is None: dic = {}
      ret = []
      n = -1
      while 1:
        n += 1
        name2 = prenam + str(n) 
        is_val, sub = subdic(dic, str(n))
        if is_val is None and n >= spyderform.length: break
        mdefault = None
        if default is not None and len(default) > n: mdefault = default[n]        
        if is_val is None:
          #absolutely no CGI entries for this index... must be non-existing in the form
          v = None
        else:  
          sf = spyderform[n] if n < spyderform.length else spyderform[0]        
          v = _reduce_cgi(sub, sf, name2, is_val, sf._typetree, mdefault)
        ret.append(v)
      pos = len(ret) - 1
      while pos >= 0 and ret[pos] == None: pos -= 1
      ret = ret[:pos+1]
      if ret == [None] * len(ret): ret = None
      v = ret            
    else:
      assert is_value is not False, (prename, spyderform.typename)      
      v = dic
      if v is not None and issubclass(cls, getattr(Spyder, "String")): v = v.replace("\r\n", "\n")
      formdefault = default
      
      if default is not None and tostr(v) == tostr(default): 
        v = None      
      elif typ == "file": #no resource
        formdefault = None #a filled-in file is always a form modification
        filetyp = getattr(Spyder, spyderform.file)
        if not hasattr(spyderform, "filename") or spyderform.filename is None:
          raise TypeError("Cannot parse CGI for %s (%s): spyderform has no attribute 'filename'" % (spyderform.typename, prenam[:-1]))
        if v is not None and len(v):
          v = {
            "name": spyderform.filename,
            "fileformat": filetyp.typename(),
            "format": filetyp.typename(),
            "mode": "w",
            "data": v,
          }          
        else:
          v = None
      else:        
        if typ in ("checkbox", "switch"):
          if v is not None: 
            if formdefault == True:
              v = None
            else: #implicit False
              formdefault = False
              v = True
          else: 
            if formdefault == True: #explicit False
              v = False
            else: #implicit False
              formdefault = False
              v = None
        else:        
          if v is None or v == "": 
            if formdefault is None or str(formdefault) == "": #implicit empty string
              v = None
            else: #explicit empty string, if we got one
              pass
  else: #composite Spyder object 
    assert is_value is not True, (prename, spyderform.typename)
    if dic is None: dic = {}
    v = {}
    for membername in spyderform.get_membernames():
      name2 = prenam + membername
      member = getattr(spyderform, membername)
      is_val, sub = subdic(dic, membername)
      mtt = [e[1] for e in typetree.members if e[0] == membername]
      if len(mtt) == 0: continue #"virtual" form member
      if len(mtt) > 1: raise Exception(membername)
      mtt = mtt[0]
      mdefault = None
      if default is not None: mdefault = getattr(default, membername)
      vv = _reduce_cgi(sub, member, name2, is_val, mtt, mdefault)
      if vv is not None: v[membername] = vv
    for k in dic:
      if k.startswith("_"): v[k] = dic[k]
    if not len(v): v = None
  return v

def dict_from_fieldstorage(fs):  
  dic = {}
  for k in fs:
    dic[k] = fs.getfirst(k)
  return dic
  
def reduce_cgi(dic, form):
  """
  Parses a CGI dict (FieldStorage object) or a plain Python dict
  The second argument can be a Spyder class, object or form

  Returns a reduced CGI dict that only contains the values different from formdefault (delta)
  """
  if isinstance(dic, stdcgi.FieldStorage):
    dic = dict_from_fieldstorage(dic)
  return _reduce_cgi(dic, form, None, False, form._typetree, None)

def filter_delta(dic, form):
  """
  Filters a delta from empty entries
  """
  return _filter_delta(dic, form, form._typetree, None)
  
def cgi(dic, spyderobj, resourceobj = None, spydertype = None):
  """
  Parses a CGI dict (FieldStorage object, or compatible object with a .getfirst() method)
  The second argument can be a Spyder class, object or form

  If a resourceobj is defined, then its resources are considered to be embedded in the model:
   as long as they haven't been redefined in the reduced dict, resources are taken from the resourcemodel
    
  Returns a Spyder object (can be None if the parsing went wrong) 
   and the model status (should be OK if the parsing went correctly)
  """
  import spyder.modules.core
  from spyder.modules.core import spyderform
  from spyder.formtools import apply_delta, model as model_, controller as controller_
  if isinstance(spyderobj, spyderform):
    form = spyderobj
  elif isinstance(spyderobj, Spyder.Object):
    form = spyderobj._form()
  else:
    ok = False
    try: 
      if issubclass(spyderobj, Spyder.Object):
        if spydertype is None: spydertype = type(spyderobj)
        form = spyderobj._form() 
        ok = True
    except TypeError:
      pass
    if not ok: raise TypeError(spyderobj)
  
  delta = reduce_cgi(dic, form)
  arraymarker = getattr(form, "arraymarker", None)
  #print("DELTA1", delta)
  if resourceobj is not None:
    update_webdict(delta, resourceobj, arraymarker)
  #print("DELTA2", delta)
  delta = filter_delta(delta, form)
  #print("DELTA3", delta)
  assert spydertype is not None
  model = model_(spydertype)
  controller = controller_(form, model)
  
  apply_delta(delta, controller)
  
  return model._value, model._status(controller), delta

def _update_webdict(webdict, m, arraymarker, parent=None, parentkey=None):
  try:
    if isinstance(m, str): raise TypeError
    items = [a for a in m]
    items = enumerate(items)
    if webdict is None:
      webdict = []
    ar = True
  except TypeError:
    items = list(m.__dict__.items())
    if webdict is None:
      webdict = {}
    ar = False
    
  #check for clonearraymarkers to support dynamic addition/deletion of arrays
  markings = {}
  if ar and len(webdict) > 0:
    for n in range(len(webdict)):
      w = webdict[n]
      if isinstance(w, dict) and arraymarker is not None and arraymarker in w:
        m = w[arraymarker]
        try:
          m = int(m)
          if m not in markings:
            markings[m] = n
        except:
          pass
  for k,v in items:
    key = k
    if len(markings) > 0:
      if k not in markings: continue
      key = markings[k]
    exists = False
    if ar:
      if len(webdict) > key: exists = True
    else:
      try:
        if k in webdict and webdict[key] is not None: exists = True
      except TypeError:
        pass
    if isinstance(v, Spyder.Resource): 
      if k in webdict and webdict[key] is not None: continue
      webdict[key] = v.dict()
    elif isinstance(v, Spyder.Object):      
      if exists:
        _update_webdict(webdict[key], v, arraymarker)
      else:
        _update_webdict(None, v, arraymarker, webdict, key)
        if len(markings) == 0:
          if ar and len(webdict) <= k: break #no more elements after the first None, except if we have clonearray markers
  
def update_webdict(webdict, m, arraymarker=None):
  """
  Updates a (expanded) CGI dict
   with resources from a different model m
  """
  return _update_webdict(webdict,m,arraymarker)  