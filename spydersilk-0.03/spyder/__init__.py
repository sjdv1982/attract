"""
front-end for loading Spyder modules
"""

import functools, imp, sys, os, traceback

spyderdir = os.path.split(__file__)[0]
if not len(spyderdir): spyderdir = "."

python2 = (sys.version_info[0] == 2)
python3 = (sys.version_info[0] == 3)

class SpyderTypeImportError(ImportError): pass
class SpyderModuleImportError(ImportError): pass

if python2:
  from .safe_eval2 import safe_eval
elif python3:  
  from .safe_eval3 import safe_eval


__spydermodules__ = {}
__current = set()
  
silent = False
corecode = None
"""
records all executed code up to the first encountered .spy file
This code is considered core code; change in this code will re-compile
all 
After the first encountered .spy file, the core is locked and the checksum computed
"""
corelock = False

def core_append(code):
  global corecode
  if corecode is None: corecode = ""
  corecode += code
  
def get_corechecksum(): return corechecksum

def core_locked(): return corelock

def lock_core():
  global corelock, corechecksum
  if corelock == False:
    import zlib
    cc = corecode
    if cc is None: cc = ""
    c = cc.encode("UTF-8")
    corechecksum = zlib.adler32(c) & 0xffffffff
    corelock = True


def __constructor(__funcname__, __class__, __defaultconstructor__, *args, **kargs):    
  try:
    __class__.__constructor__ = __funcname__
    ret = __class__(*args, **kargs)    
  finally:
    __class__.__constructor__ = __defaultconstructor__  
  return ret


def set_module(name, fullname, mod):
  if name not in globals(): 
    globals()[name] = mod
    setattr(spyder_reduced, name, mod)
  __spydermodules__[fullname] = mod

class Object(object): pass
__types__ = {"Object":Object}
arraymax = 2

from . import loader
sys.meta_path = [m for m in sys.meta_path \
 if not hasattr(m, "__name__") or not hasattr(m, "__module__") \
   or m.__name__ != "finder" or m.__module__ != loader.__name__
]  
sys.meta_path.insert(0, loader.finder)    
from .define import *

newnamespace = functools.partial(imp.new_module, "Spydernamespace")

def __exceptionhook(type,value,tb):
  import traceback
  traceback.print_exception(type,value,tb)  
  
sys.excepthook = __exceptionhook
  
errorpath = []  
def format_errorpath(errpath=None):
  if errpath is None: errpath = errorpath
  pathstr = ""
  for parentclass,memberclass,m in errpath:    
    if isinstance(m, int): pathstr += "[%d]" % m 
    elif len(pathstr) == 0: pathstr = str(m)
    else: pathstr += "."+str(m)
  if not len(pathstr): return None  
  pathstr = "\nMember %s:\n" % pathstr      
  return pathstr
  
def exception(exceptionstring=None,newline=False):  
  if exceptionstring is None:
    exc = sys.exc_info()
    if exc[0] is None: return None
    if issubclass(exc[0], ValidationError) or issubclass(exc[0], ConstructionError):    
      s = str(exc[1])
    else:
      s = traceback.format_exc()  
  else:
    s = exceptionstring    
  pathstr = format_errorpath()
  ss = s.lstrip("\n").split("\n")
  if len(ss) and not len(ss[0]): ss = ss[1:]  
  if pathstr is None: 
    ret = "\n" + "\n".join([l for l in ss])  
    if newline and (not len(ss) or len(ss[-1])): ret += "\n"
    return ret
  if len(ss) and "\n" + ss[0].strip() + "\n" == pathstr: 
    if newline and (not len(ss) or len(ss[-1])): s += "\n"
    return "\n" + s
  s2 = "\n".join(["  " + l for l in ss])
  if newline and (not len(ss) or len(ss[-1])): s2 += "\n"
  return pathstr + s2
  
from .loader import delete_wisdom, write_wisdom, export_wisdom, wisdomdict

def __load_module(module):
  from . import loader
  if not isinstance(loader.wisdom, wisdomdict): 
    loader.wisdom = wisdomdict(loader.wisdom)

  from .loader import wisdom
  if module in __current: return
  assert module in wisdom, (module, wisdom.keys())
  w = wisdom[module]
  if w.fullname in __current: return
  
  __current.add(w.fullname)
  
  if w.fullname in loader.loaded: 
    ret =__spydermodules__[w.fullname]
    if w.fullname not in sys.modules: sys.modules[w.fullname] = ret
    nam = w.fullname
    mod = sys.modules[nam]
    while nam.find(".") > -1:
      pos = nam.rindex(".")
      tail = nam[pos+1:]
      nam = nam[:pos]
      if nam not in sys.modules: break
      setattr(sys.modules[nam], tail, mod)
      mod = sys.modules[nam]
    
    __current.remove(w.fullname)
    return ret
  if not corelock:
    for m in w.coremodules: __load_module(m)
  loadtypes = set()
  loadtypes.update(w.needed) 
  loadtypes.update(w.used)
  for t in loadtypes:
    try:
      __load_type(t)
    except SpyderTypeImportError as e:
      raise SpyderModuleImportError(module, e.args[0])
  old_tempdir = loader.tempdir  
  old_recompile = loader.recompile
  loader.tempdir, loader.recompile = w.tempdir, w.recompile
  l = loader.loader(w.filename, w.path, w.name)
  ret = l.load_module(w.fullname,check_local=False)
  loader.tempdir, loader.recompile = old_tempdir, old_recompile
  __current.remove(w.fullname)
  return ret

typedep = {}
def __load_type(t):
  if t == "Object" or t == "PathList": return
  if t in __current: return
  parse_wisdom()
  __current.add(t)
  t0 = t
  while t0.endswith("Array"): t0 = t0[:-len("Array")]
  if t0 != "Resource" and t0.startswith("Resource"): t0 = t0[len("Resource"):]
  if t0 not in typedep:
    raise SpyderTypeImportError(t)
  __load_module(typedep[t0])
  __current.remove(t)
  return 

parsed_wisdom = False
def parse_wisdom():
  global parsed_wisdom
  from .loader import wisdom  
  if not parsed_wisdom:
    wis = sorted(list(wisdom.values()),key=lambda w: w.priority)
    for w in wis:
      for t in w.defined: typedep[t] = w.name
    parsed_wisdom = True

def load(obj):
  parse_wisdom()
  from .loader import wisdom    
  if obj == "*":
    ret = []
    for w in wisdom.keys(): ret.append(__load_module(w))
  elif obj == "__all__": #all types
    ret = {}
    for t in typedep: 
      v = __load_type(t)
    ret = __types__
  elif spydercompile.validvar2(obj) or obj == "Object":
    if obj not in __types__:
      __load_type(obj)
    ret = __types__[obj]
  else:
    ret = __load_module(obj)
  return ret    

def _relative_import(modulename):
  oldmod = None
  if modulename in sys.modules: 
    oldmod = sys.modules[modulename]
    del sys.modules[modulename]
  oldpath = list(sys.path)
  sys.path.append(os.getcwd())
  try:
    file, pathname, desc = imp.find_module(modulename, [os.getcwd()])
    mod = imp.load_module(modulename, file, pathname, desc)
  finally:
    sys.path = oldpath
    if oldmod is not None:
      sys.modules[modulename] = oldmod
    elif modulename in sys.modules:
      del sys.modules[modulename]
  return mod

from .spydercompile import validvar, validvar2

from .moduleSpyder_wrapper import moduleSpyder_wrapper
sys.modules["Spyder"] = moduleSpyder_wrapper()

spyder_reduced = imp.new_module("spyder")
forbidden = ["functools", "imp", "sys", "os", "traceback", 
 "typedep", "core_append", "lock_core", 
 "set_module", "loader", "delete_wisdom", "write_wisdom", "export_wisdom", "wisdomdict"]
 
class ConstructionError(Exception): pass
class ValidationError(Exception): pass

for k in list(globals().keys()):
  if k not in forbidden: setattr(spyder_reduced,k, globals()[k])
spyder_reduced._file_load =  loader.file_load
spyder_reduced.linesep = os.linesep

if len(loader.wisdom) == 0 and not '--make-wisdom' in sys.argv:
  #We are in a developer environment, no pre-generated WISDOM filename
  import spyder.modules.core
  import spyder.modules.builtin
  import spyder.modules.basic
  import spyder.modules.atom
  import spyder.modules.http
  import spyder.modules.tarantula
  import spyder.modules.models3d
  import spyder.modules.canvas

