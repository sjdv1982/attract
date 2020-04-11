from __future__ import print_function

import imp, sys,os, tempfile,functools

"""
The following functions are used to retrieve local files.
They are used by this loader (.spy and .fly files, wisdom files)
The _file_... functions are used as well as by core.parse (fromfile) and builtin.File
Redirect them if you want to implement additional file loading mechanisms 
"""
_file_exists = os.path.exists
_file_load = open
_change_dir = os.chdir
_file_access = os.access

def file_exists(filename): 
  return _file_exists(filename)

def file_load(filename, mode = None):
  if mode is None: return _file_load(filename)
  return _file_load(filename, mode)

def file_access(filename, okmods):
  return _file_access(filename, okmods)
    
def change_dir(directory): 
  return _change_dir(directory)
  
from . import __types__, python2, python3
from . import lock_core, core_locked, core_append, get_corechecksum
from .spydercompile import spydercompile, update_dependencies, cache_from_source

if python2:
  try:
    import cPickle as pickle
  except ImportError:
    import pickle
else:
  import pickle    
    
tempdir = tempfile.gettempdir()

#recompile can have three values:
# - True: Always recompile
# - False: Never recompile
# - None: Recompile if changed
recompile = None 

typedefs = {}
coremodules = []
loaded = set()
relpathmode = False

class Wisdom(object):
  def __init__(self, filename, path, name, fullname, defined, needed, used):
    priority = 1
    if len(wisdom.keys()):
      priority = max([w.priority for w in wisdom.values()]) + 1
    self.priority = priority
    self.coremodules = list(coremodules)
    
    self.filename = os.path.abspath(filename)
    self.path = path
    self.name = name
    self.fullname = fullname
    self.defined = set(defined)
    self.needed = set(needed)
    self.used = set(used)
    self.tempdir = tempdir
    self.recompile = recompile
    self.is_local = False

wisdomfile = os.path.split(__file__)[0] + os.sep + "WISDOM"

def substitute_path(p, path):
  pp = p.replace("#PATH", path)
  if not pp.startswith("#"): pp  = os.path.abspath(pp)
  return pp

def substitute_spyderdir(p):
  from . import spyderdir
  pp = p.replace("#SPYDERDIR", spyderdir)
  if not pp.startswith("#"): pp  = os.path.abspath(pp)
  return pp
  
def desubstitute_spyderdir(p):
  from . import spyderdir
  p = substitute_spyderdir(p)
  if p.startswith(spyderdir): 
    p = "#SPYDERDIR" + p[len(spyderdir):]
  elif relpathmode == True:
    p = "#SPYDERDIR" + os.sep + os.path.relpath(p, spyderdir)
  return p  
  

class wisdomdict(dict):
  def __getitem__(self, key):
    try:
      return dict.__getitem__(self, key)
    except KeyError:
      ret = []
      for w in sorted(self.values(),key=lambda x:str(x.name)):
        if w.name == key: ret.append(w)
      if not len(ret): raise
      assert len(ret) == 1, (len(ret), key)
      return ret[0]
  def __contains__(self, key):
    if dict.__contains__(self, key): return True
    for w in self.values():
      if w.name == key: return True
    return False

wisdom = wisdomdict()
try:
  w = file_load(wisdomfile, "rb")
  wisdom = pickle.load(w)
  if not isinstance(wisdom, wisdomdict):
    wisdom = wisdomdict(wisdom)
except:
  pass
for wis in wisdom.values():
  if wis.path is not None:
    for pnr, p in enumerate(wis.path):
      wis.path[pnr] = substitute_spyderdir(p)
  wis.filename = substitute_spyderdir(wis.filename)

def default_core_wisdom():
  from . import load, corecode
  global wisdom
  if not isinstance(wisdom, wisdomdict):
    wisdom = wisdomdict(wisdom)
  if corecode is not None: 
    wisdom["default-core"] = Wisdom("",None,None,None,[],[],[])
    return
  if "default-core" in wisdom:
    for m in wisdom["default-core"].coremodules:
      load(m)
  

def delete_wisdom(key=None):
  global wisdom
  if key is None:
    wisdom = wisdomdict()
  else:
    wisdom.pop(key, None)  

def load_needed(filename, needed):
  from . import load, parse_wisdom, typedep
  parse_wisdom()
  for typ in sorted(needed):
    if typ not in typedep:
      raise ImportError("'%s' needs unknown Spyder type '%s'" % (filename, typ))
    load(typ)
  
def write_wisdom():
  default_core_wisdom()
  global spyderdir
  from . import spyderdir
  try:
    w = file_load(wisdomfile, "wb")
  except:
    raise Exception("Unable to open wisdom file for writing; you may need root/admin privileges to write wisdom")
  for wis in wisdom.values():
    if wis.is_local: continue
    if wis.path is not None:
      for pnr, p in enumerate(wis.path):
        wis.path[pnr] = desubstitute_spyderdir(p)
    wis.filename = desubstitute_spyderdir(wis.filename)
  pickle.dump(wisdom, w, pickle.HIGHEST_PROTOCOL)  
  w.close()

def export_wisdom(key):
  default_core_wisdom()
  wis = wisdom[key]
  global spyderdir
  from . import spyderdir
  wis.filename = desubstitute_spyderdir(wis.filename)
  localwisdomfile = os.path.split(wis.filename)[0] + os.sep + "WISDOM"
  try:
    w = file_load(localwisdomfile, "wb")
  except:
    raise Exception("Unable to open local wisdom file '%s' for writing; you may need root/admin privileges to write wisdom" % localwisdomfile)
  if wis.path is not None:
    for pnr, p in enumerate(wis.path):
      wis.path[pnr] = desubstitute_spyderdir(p)
  pickle.dump(wis, w, pickle.HIGHEST_PROTOCOL)  
  w.close()

def _load_python(f,modulename):
  from .spydercompile import compile_obj, get_types
  ff = file_load(f)
  code = ff.read()
  ff.close()
  bytecode = None
  f2 = cache_from_source(f)
  if file_exists(f2):  
    ff2 = file_load(f2,"rb")
    bytecode = ff2.read()
    ff2.close()
  obj, bytecode_ok = compile_obj(code,modulename+"::"+f,f,bytecode,recompile)
  typef = f[:-len(".py")] + ".TYPES"
  defined, dmmy, used = get_types(f, code, typef, True, cache=True)
  for t in defined:
    if t not in typedefs: typedefs[t] = []  
  if not bytecode_ok:
    if bytecode is not None:
      try:
        os.remove(f2)
      except:
        pass
    import py_compile
    try:
      py_compile.compile(f,cfile=f2,dfile=modulename+"::"+f)
    except:
      pass      
  return code, obj, defined, used

header = """
from __future__ import print_function
from spyder import spyder_reduced as spyder
import Spyder, functools
"""

class finder(object):
  @staticmethod
  def find_module(fullname,path=None):
    try:
      os.sep
    except AttributeError: #Python clean up phase, just leave...
      return  
    name = fullname.split(".")[-1]
    paths = path
    if path is None: paths = sys.path
    for p in paths:
      if p == "": p = "."
      fnams =(p+os.sep+name+".fly",
              p+os.sep+name+".spy",
              p+os.sep+name+os.sep+"__init__.fly"
             )
      for fnam in fnams:       
        if file_exists(fnam): 
          return loader(fnam,path,name)

def flyloader(flyname, modulename, path, temp_prefix):
  definedtypes, neededtypes, usedtypes = set(), set(), set()

  code = []
  fil = file_load(flyname)
  olddir = os.getcwd()
  change_dir(os.path.split(flyname)[0])
  flines = fil.read().splitlines()
  for linenr,f in enumerate(flines):
    f = f.replace('\r\n', "")
    f = f.replace('\n', "")
    f = f.strip()
    if len(f) == 0: continue
    if f.startswith("#"): continue
    if not file_exists(f):
      raise ImportError("Error in %s, line %d: file %s does not exist" % \
       (flyname, linenr, f))
    if f.endswith(".spy"):
      fdir = os.path.split(os.path.abspath(f))[0]
      default_core_wisdom()
      lock_core()
      s = file_load(f).read()
      corechecksum = get_corechecksum()
      tempdir2 = substitute_path(tempdir, path)
      c, obj, defined, needed, used = spydercompile(
        s,tempdir2,recompile,
        temp_prefix,corechecksum,typedefs,modulename,f
      )
      definedtypes.update(defined)     
      if needed is not None: neededtypes.update(needed)
      usedtypes.update(used)
      code.append((fdir,f,c,obj))     
    elif f.endswith(".fly"):
      flynam = os.path.split(f)[1][:-len(".fly")]
      new_temp_prefix = temp_prefix + "_" + flynam
      newcode, defined, needed, used = flyloader(f, modulename, path, new_temp_prefix)
      definedtypes.update(defined)
      neededtypes.update(needed)
      usedtypes.update(used)
      code += newcode
    elif f.endswith(".py"):
      fdir = os.path.split(os.path.abspath(f))[0]
      c, obj, defined, used = _load_python(f,modulename)
      definedtypes.update(defined)
      usedtypes.update(used)
      code.append((fdir, f, c,obj))
  fil.close()
  change_dir(olddir)
  update_dependencies(definedtypes, neededtypes, usedtypes)
  return code, definedtypes, neededtypes, usedtypes

class loader(object):
  def __init__(self, filename, path, name):
    self.filename = filename
    self.path = path
    self.name = name
    self.moddicts = {}
      
  def load_module_singlefile(self, filedir, fullname, code, obj):
    mod = sys.modules.setdefault(fullname, imp.new_module(fullname))
    mod.__dict__.update(__types__)
    mod.__file__ = self.filename
    mod.__loader__ = self
    if self.path is not None:
      mod.__path__ = self.path

    olddir = os.getcwd()
    change_dir(filedir)
    ok = True
    self.currcode = code
    exec(header, mod.__dict__)
    exec(obj, mod.__dict__)
    change_dir(olddir)

    if not core_locked(): core_append(code)
    from . import set_module
    set_module(self.name,fullname,mod)
    return mod  
  
  def load_module(self,fullname,check_local=True):        
    global wisdom
    if not isinstance(wisdom, wisdomdict):
      wisdom = wisdomdict(wisdom)
    
    filedir = os.path.split(os.path.abspath(self.filename))[0]    
    localwisdomfile = filedir + os.sep + "WISDOM"
    if check_local and file_exists(localwisdomfile):
      w = file_load(localwisdomfile, "rb")
      wis = pickle.load(w)
      assert isinstance(wis, Wisdom), localwisdomfile
      wis.path = self.path
      if wis.path is not None:
        for pnr, p in enumerate(wis.path):
          wis.path[pnr] = substitute_spyderdir(p)
      wis.filename = substitute_spyderdir(self.filename)
      wisdom[fullname] = wis 
      from . import load
      return load(fullname)

    #Load a single file
    if self.filename.endswith(".spy"):
      default_core_wisdom()
      lock_core()
      s = file_load(self.filename).read()
      corechecksum = get_corechecksum()
      tempdir2 = substitute_path(tempdir, filedir)
      code, obj, defined, needed, used = spydercompile(
       s,tempdir2,recompile,
       fullname,corechecksum,typedefs,self.name,self.filename
      )
      load_needed(self.filename, needed)
      if not core_locked(): coremodules.append(self.name)
      wisdom[fullname] = Wisdom(self.filename, self.path, self.name, fullname, defined, needed, used)
      from . import typedep
      for t in defined: typedep[t] = self.name

      ret = self.load_module_singlefile(filedir, fullname, code, obj)
      loaded.add(fullname)
      return ret

    #Load a file list          
    allcode = ""
    code, defined, needed, used = flyloader(self.filename, fullname,filedir,self.name)
    load_needed(self.filename, needed)
    if not core_locked(): coremodules.append(self.name)
    wisdom[fullname] = Wisdom(self.filename, self.path, self.name, fullname, defined, needed, used)
    from . import typedep
    for t in defined: typedep[t] = self.name
    
    mod = sys.modules.setdefault(fullname, imp.new_module(fullname))    
    mod.__file__ = self.filename
    mod.__loader__ = self
    if self.path is not None:
      mod.__path__ = self.path

    for fdir, fnam, currcode,obj in code:
      self.currcode = currcode
      olddir = os.getcwd()
      change_dir(fdir)
      ok = True
      mod.__dict__.update(__types__)    
      exec(header, mod.__dict__)
      mdict = dict(mod.__dict__)
      mdict["__loader__"] = fakeloader(fnam, self.currcode)
      exec(obj, mdict)
      del mdict["__loader__"]
      mod.__dict__.update(mdict)
      for fnam0 in self.moddicts:
        self.moddicts[fnam0].update(mdict)
      mdict["__loader__"] = fakeloader(fnam, self.currcode)
      self.moddicts[fnam] = mdict
      mod.__loader__ = self
      change_dir(olddir)
      allcode += currcode
    self.currcode = allcode
    if not core_locked(): core_append(allcode)
    from . import set_module
    set_module(self.name,fullname,mod)
    loaded.add(fullname)
    return mod  
  def get_source(self,fullname):
    return self.currcode
    
class fakeloader(object):
  def __init__(self, name, code):
    self.name = name
    if hasattr(code, "decode"):
      code = code.decode("UTF-8")
    self.code = code
  def get_source(self,fullname):
    return self.code
    
