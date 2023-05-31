import copy, spyder
spyderforms = {}

maxlength_default = 1000

class dummyclass(object):
  def __setattr__(self, attr, value):
    pass
  def __getattr__(self, attr):
    return dummyclass()
  def __call__(self, *args):
    return dummyclass()

class groupelement(object):
  def __init__(self, id, type):
    self.id = id
    self.type = type
    self.members = []
    

class spyderform_othertokenwrapper(object):
  def __init__(self, attr, othertokenlist):
    self._attr = attr
    self._othertokenlist = othertokenlist
  def __call__(self, *args):
    self._othertokenlist.append((self._attr, args))

class spyderform(object):
  typetree_attrs = ("typename", "arraycount", "is_resource", "is_default")
  _members = None
  _membernames = None
  _arrayforms = None  
  _props = {}
  def __init__(self, typetree = None, toplevel=True, element=False):
    if typetree is None: return
    nam = typetree.typename
    if typetree.is_resource: nam = "Resource" + nam
    self._typetree = typetree    
    self._props = {}
    self._props["is_resource"] = typetree.is_resource
    self._props["default"] = typetree.default    
    self._props["default_expr"] = typetree.default_expr
    for attr in self.typetree_attrs:
      self._props[attr] = getattr(typetree, attr)
    if hasattr(typetree, "spydertype"): 
      self._props["spydertype"] = typetree.spydertype
    if hasattr(typetree, "type") and typetree.type is not None: 
      if not hasattr(typetree, "typemapped") or typetree.typemapped == False:
        self._props["type"] = typetree.type
      else:
        self._props["type"] = typetree.typename
        
    if toplevel: 
      del self._props["default"]
      del self._props["is_default"]
    self._groups = []
    if element: self._props["arraycount"] -= 1
    if self._props["arraycount"] > 0: 

      self._props["maxlength"] = maxlength_default
      if hasattr(typetree, "maxformlength"): 
        self._props["maxlength"] = typetree.maxlength

      self._arrayforms = {}
      typ = typetree.typename+"Array"*(self._props["arraycount"]-1)
      if self._props["is_resource"]: typ = "Resource" + typ
      if typ not in spyderforms:
        typ2 = typetree.typename+"Array"*self._props["arraycount"]
        if self._props["is_resource"]: typ2 = "Resource" + typ2
        assert typ2 in spyder.__types__, "Unknown Spyder type %s" % typ2
        typtree = spyder.__types__[typ2]._typetree()
        self._ele = spyderform(typtree, toplevel=False, element=True)
      else:
        form2 = spyderforms[typ]      
        self._ele = form2.get_copy()
      self._ele._props["default"] = None
      self._ele._props["is_default"] = False

    if element: return
    if self._props["arraycount"] == 0 and typetree.members is not None:
      self._members = {}
      self._membernames = [m[0] for m in typetree.members] 
      for name, membertree in typetree.members:
        memtype = None
        if membertree.typename is not None:
          memtype = membertree.typename + "Array" * membertree.arraycount
        if membertree.is_resource: memtype = "Resource" + memtype
        if membertree.members is None and memtype not in spyderforms:
          spyderform(membertree)

        if membertree.arraycount and memtype not in spyderforms:
          spyderform(membertree)  

        if memtype is None:
          v = spyderform(membertree)
        else:  
          v = spyderforms[memtype].get_copy()
        v._props["default"] = membertree.default
        for attr in self.typetree_attrs:
          v._props[attr] = getattr(membertree, attr)        
        self._members[name] = v
    if toplevel and typetree.typename is not None: 
      name = typetree.typename+"Array"*typetree.arraycount
      if typetree.is_resource: name = "Resource" + name
      spyderforms[name] = self
  
  def new_group(self, id, type):
    g = groupelement(id, type)
    self._groups.append(g)
    return g      
  def get_arrayforms(self):
    assert self._props["arraycount"] > 0
    return self._arrayforms.copy()
  def get_properties(self):
    return self._props.copy()
  def get_membernames(self):
    if self._membernames is None: return None
    ret = []
    if "memberorder" in self._props and self._props["memberorder"] is not None:
      for m in self._props["memberorder"]:
        if m in self._membernames: ret.append(m)
      for m in self._membernames:
        if m not in ret: ret.append(m)        
      return ret
    else:
      return list(self._membernames)
  def get_copy(self):
    ret = spyderform()
    ret._typetree = self._typetree
    ret._props = self._props.copy()
    ret._arrayforms = None
    if self._arrayforms is not None:
      ret._arrayforms = dict([(key, value.get_copy()) for key,value in self._arrayforms.items()])
      ret._ele = self._ele.get_copy()
    ret._members = None
    if self._members is not None:
      ret._members = dict([(key, value.get_copy()) for key,value in self._members.items()])
    ret._membernames = self._membernames      
    ret._groups = list(self._groups)
    return ret      
  def __setattr__(self, attr, value):
    if attr in ("default", "value"):
      if not isinstance(value, tuple):
        value = (value,)
      if hasattr(self._typetree, "spydertype"): 
        cls = self._typetree.spydertype
      else: 
        name = self.typename+self.arraycount*"Array"
        if self._props["is_resource"]: name = "Resource" + name
        cls = spyder.__types__[name]
      if len(value) == 1 and (isinstance(value[0], cls) or value[0] is None):
        value = value[0]
      else:
        value = cls(*value)
      self._props[attr] = value  
      self._props[attr+"_expr"] = repr(value)
      return
      
    if attr.startswith("_"): 
      self.__dict__[attr] = value
      return
    
    if self._members is not None and attr in self._members:
      raise TypeError("Cannot set form property %s of %s since it is a data member" % (attr, self._typetree.typename))
    elif attr in self.typetree_attrs:
      raise TypeError("Cannot set form property %s of %s since it is constant" % (attr, self._typetree.typename))
    self._props[attr] = value
  def __getattr__(self, attr):
    if self._members is not None and attr in self._members: 
      return self._members[attr]
    elif attr.startswith("add_"):
      attr = attr[len("add_"):]
      if "_othertokens" not in self._props: self._props["_othertokens"] = []    
      return spyderform_othertokenwrapper(attr, self._props["_othertokens"])
    elif attr.startswith("get_"):  
      attr = attr[len("get_"):]
      def getter():
        if "_othertokens" not in self._props: return []
        match = [v[1] for v in self._props["_othertokens"] if v[0] == attr]
        ret = tuple()
        for m in match: ret = ret + m
        return ret
      return getter  
    elif attr in self._props: 
      return self._props[attr]
    else:
      raise AttributeError

  def __getitem__(self, key):
    if self._arrayforms is None: raise TypeError
    if key is None: return self._ele
    if not isinstance(key, int): raise TypeError(key)
    if "length" in self._props:
      length = self._props["length"]
      if key < 0: key = length + key #length - -key
      if key >= length: raise IndexError(key)
    maxlen = self._props["maxlength"]
    if key >= maxlen or key < -maxlen: raise IndexError(key)
    if key not in self._arrayforms:
      v = self._ele.get_copy()
      if "name" not in v._props:
        k = key
        if "count_from_one" in self._props \
         and self._props["count_from_one"]:
           k += 1
        if "name" in self._props:
          v._props["name"] = self._props["name"] + " " + str(k)
        else:
          v._props["name"] = str(k)
      self._arrayforms[key] = v
    return self._arrayforms[key]
  def __contains__(self, key):
    return key in self._arrayforms
  def __delitem__(self, key):
    if key in self._arrayforms: del self._arrayforms[key]
  def __iter__(self):
    ret1 = sorted([k for k in self._arrayforms.keys() if k >= 0])
    ret2 = list(reversed(sorted([k for k in self._arrayforms.keys() if k < 0])))
    return iter(ret1 + ret2)
  def __str__(self):
    if self._arrayforms is not None: 
      ret = [self._arrayforms[k] for k in self.__iter__]
      ret = [str(v) if isinstance(v, spyderform) else v for v in ret]
    else:
      ret = self._members
      if ret is not None:
        ret = ret.copy()
        for k,v in ret.items():
          if isinstance(v, spyderform):
            ret[k] = str(v)
    return str(ret)
def generate_form(typename, parentnames, source, members, deleted_members, block):
  requiredmembers, defaultmembers, optionalmembers, args, allargs = spyder.core.parse_members(typename,members,None, spyder.safe_eval)  
  s = """  
  @classmethod
  def _form(cls):
    \"\"\"Auto-generated by Spyder:
     module core
     file form.py
     function generate_form
     Returns the spyderform of the class
     This function is automatically called at class creation\"\"\"
    return spyder.core.spyderforms[cls.typename()]
  @classmethod
  def __form__(cls):
    \"\"\"Auto-generated by Spyder:
     module core
     file form.py
     function generate_form
     Generates the spyderform of the class
     This function is automatically called at class creation\"\"\"
"""
  ss = "    self = spyder.core.spyderform(cls._typetree())\n"    
  if block is not None:
    for m in deleted_members:
      ss += "    %s = spyder.core.dummyclass()\n" % (m)
    for m in requiredmembers+defaultmembers:
      ss += "    %s = self.%s\n" % (m[1],m[1]) 

    spaces = -1
    for l in block.split('\n'):
      if len(l.strip()) == 0: continue
      if spaces == -1:
        spaces = len(l) - len(l.lstrip())
        s += ss
      s += "    %s\n" % l[spaces:]    
    if spaces == -1: s += "    pass\n"
  else:
    s += ss
  s2 = """  
  @classmethod
  def _form(cls):
    \"\"\"Auto-generated by Spyder:
     module core
     file form.py
     function generate_form
     Returns the spyderform of the class
     This function is automatically called at class creation\"\"\"
    return spyder.core.spyderforms[cls.typename()]
  @classmethod
  def __form__(cls):
    \"\"\"Auto-generated by Spyder:
     module core
     file form.py
     function generate_form
     Generates the spyderform of the class
     This function is automatically called at class creation\"\"\"
    spyder.core.spyderform(cls._typetree()) 
    t = cls.typename()
    while t.endswith("Array"): t = t[:-len("Array")]
    spyder.core.spyderforms[cls.typename()]._members = spyder.core.spyderforms[t]._members
"""

  return s,s2

spyder.defineunimethod("form", generate_form)

