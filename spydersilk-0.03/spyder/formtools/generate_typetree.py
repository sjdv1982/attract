"""
generate:
Generates a typetree out of arbitrary tuples and dicts
tuples must be (v1,v2,v3,...) 
dicts must be {name1:v1,name2:v2,name3:v3,...}

names must be strings
v must be: 
- float,int,str,bool
- or: "float","str","int","bool"
- or: a Spyder type or typename
- or: another tuple or dict (nested input is supported)
"""

import copy
import Spyder, spyder
from Spyder import Object
from spyder import validvar2

from ..spydercompile import reserved_membernames

def _generate_typetree(members,typemap={}):
  #adapted from spyder.core.get_typetreedict

  ret = spyder.core.typetreeclass()
  ret.typename = None
  ret.type = None
  ret.arraycount = 0
  ret.is_default = False
  ret.is_resource = False
  ret.default = None
  ret.default_expr = None
  ret.members = []
  ret.membernames = [m[1] for m in members]
  
  for mem in members:
    is_default, mname, mtypename, mtype, default = mem[:5]
    fulltypename = None
    if len(mem) == 6: fulltypename = mem[5]
    if isinstance(mtype, dict):
      mtree = generate_typetree(mtype)
    elif not validvar2(mtypename):
      mtree = copy.copy(spyder.core.get_typetreedict("String"))
      mtree.typename = mtypename
      mtree.type = typemap[mtypename]
      mtree.typemapped = True
      mtree.spydertype = mtree.type
    else:
      mtree = copy.copy(spyder.core.get_typetreedict(mtypename))  
    mtree.is_default = is_default
    if fulltypename is not None:
      mtree.fulltypename = fulltypename
    if default is not None:
     mtree.default_expr = default
     if not validvar2(mtypename):
       try:
         mtree.default = mtree.type(default)
       except:
         print("Expression: %s" % str(default))
         raise
     elif mtypename in ("String", str):
       mtree.default_expr = default
     elif not isinstance(default, float) \
      and not isinstance(default, int) \
      and not isinstance(default, bool):
       try:
         mtree.default_expr = eval(str(default), spyder.__types__.copy())
       except:
         print("Expression: %s" % str(default))
         raise
     if mtree.default_expr is not None:
       var = "__default_expr__" 
       dic = spyder.__types__.copy()
       dic[var] = mtree.default_expr
       exp = "%s(%s)" % (mtypename, mtree.default_expr)       
       try:
         mtree.default = eval("%s(%s)" % (mtypename, var), dic)
       except Exception as e:
         print("Expression: %s" % exp)
         raise
    ret.members.append((mname, mtree))  
  
  return ret
   
spydermap = {
  "str" : "String",
  "bool" : "Bool",
  "int" : "Integer",
  "float" : "Float", 
  "number": "Float",
  "real": "Float",
}  
  
def generate_typetree(template, typemap={}):
  if isinstance(template, list):
    for it in template:
      assert isinstance(it, tuple), it
      assert len(it) == 2, it
    members = tuple(template)  
  elif isinstance(template, tuple):
    members = tuple([("par%d" % (n+1),t) for n,t in enumerate(template)])
  elif isinstance(template, dict):
    for k in template: assert isinstance(k, str), k
    members = list(template.items())
    members.sort(key=lambda i: i[0])
    members = tuple(members)
  else:
    raise TypeError("input must be tuple or dict")
  
  assert len(members)
  
  members2 = []  
  for membername,t in members:    
    #in case of tuples, we're assuming hive system style annotated type names;
    #- nested templates are not supported!
    #- if t is a 2-tuple, the 2nd item is interpreted as default value
    
    assert membername not in reserved_membernames, membername
    assert not membername.startswith("_"), membername
    default = None
    if isinstance(t, tuple) and len(t) == 2:
      t, default = t
    fulltypename = None
    if isinstance(t, tuple): 
      if all((isinstance(v, str) for v in t)):
        fulltypename = t
    while isinstance(t, tuple): t = t[0] 
    
    if isinstance(t, str):
      mtypename = t
      if validvar2(t):
        mtype = getattr(Spyder, t)
      elif t in ("str","bool","int","float"):
        mtypename = spydermap[t]
        mtype = getattr(Spyder, mtypename)
      elif t in typemap:
        mtype = typemap[t]
      else: 
        raise TypeError("Unknown type", t)
    elif isinstance(t, list):    
      raise TypeError("Generating typetrees from lists is not supported", t)
    elif isinstance(t, dict):  
      mtype = t
      mtypename = None
    else:
      assert callable(t), t
      mtype = t
      if issubclass(t, Object):
        mtypename = t.typename()
      elif t in (str, bool, int, float):
        mtypename = spydermap[t.__name__]
        mtype = getattr(Spyder, mtypename)
      else: 
        mtypename = None
    is_default = default is not None    
    mem = (is_default, membername, mtypename, mtype, default)
    if fulltypename is not None:
      mem = (is_default, membername, mtypename, mtype, default, fulltypename)
    members2.append(mem)     
  
  tree = _generate_typetree(members2,typemap)
  return tree
  


