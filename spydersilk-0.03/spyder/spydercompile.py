from __future__ import print_function

import os, zlib, ast, imp, marshal

import re, sys
from . import safe_eval
from copy import copy
import struct

reservedtypes = ("Spyder", "Type", "Object", "Delete", "Include", "Exception", "None", "True", "False")
reserved_membernames = (
  "type","typename","subtype","convert","cast","validate",
  "data","str","repr","dict","fromfile",
 "tofile","totempfile","threadconvert",
 "listen","block","unblock","buttons","form","length","name","type"
)

def cache_from_source(f):
  try:
    ret = imp.cache_from_source(f)
    if f.endswith(".spy.py"):
      rr = ret.split(".")
      if len(rr) < 4 or rr[-3] != "spy":
        ff = os.path.split(f)[1][:-3]
        ret0 = os.path.split(ret)[0]        
        ret = ret0 + os.sep + ff + "." + rr[-2] + "." + rr[-1]        
  except:
    ret = f+"c"
  return ret  
  

def validvar(s):
  if s.replace("_", "x").isalnum() == False: return False
  if s[0].isupper() == False: return False
  if len(s) > 1 and s == s.upper(): return False
  if s.endswith("Array") or s.endswith("Error") or s.endswith("Exit") or s in reservedtypes: return False
  return True

def validvar2(s):
  if not isinstance(s, str): return False
  while s.endswith("Array"): 
    s = s[:-len("Array")]
  return validvar(s)

quotematch = re.compile(r'(([\"\']).*?\2)')
triplequotematch = re.compile(r'(\"\"\"[\w\Wn]*?\"\"\")')
curlymatch = re.compile(r'{[^{}]*?}')

def getlines(s):
  """
  Parses a single logic line from a Spyder context 's' 
  Within a file context, every Type, Define, etc. is a single line
  Within a Type definition context, every every block definition and member definition is a single line
  Triple-quoted strings are automatically removed
  """
  pos = 0
  s0 = ""
  for pp in triplequotematch.finditer(s):
    s0 += s[pos:pp.start()] + "&" * (pp.end()-pp.start())
    pos = pp.end()
  s0 += s[pos:]

  p = quotematch.finditer(s0)
  pos = 0
  mask0 = ""
  for pp in p:
      mask0 += s[pos:pp.start()] + "*" * (pp.end()-pp.start())
      pos = pp.end()
  mask0 += s[pos:]

  while 1:
    p = curlymatch.finditer(mask0)
    pos = 0
    mask = ""
    for pp in p:
      mask += mask0[pos:pp.start()] + "!" * (pp.end()-pp.start())
      pos = pp.end()
    mask += mask0[pos:]
    if pos == 0: break
    mask0 = mask

  
  pos = 0
  mask = ""
  for pp in triplequotematch.finditer(mask0):
    mask += mask0[pos:pp.start()] + "!" * (pp.end()-pp.start())
    pos = pp.end()
  mask += mask0[pos:]
  
  lines0 = mask.split("\n")
  lines = []
  pos = 0
  for l in lines0:
    pos2 = pos + len(l)
    lines.append(s[pos:pos2])
    pos = pos2 + len("\n")
  return lines

def getblock(s):
  
  pos = 0
  s0 = ""
  for pp in triplequotematch.finditer(s):
    s0 += s[pos:pp.start()] + "&" * (pp.end()-pp.start())
    pos = pp.end()
  s0 += s[pos:]
  
  p = quotematch.finditer(s0)
  pos = 0
  mask0 = ""
  for pp in p:
      mask0 += s[pos:pp.start()] + "*" * (pp.end()-pp.start())
      pos = pp.end()
  mask0 += s[pos:]

  preblock = s
  postblock = ""
  blocks = []
  while 1:
    p = curlymatch.finditer(mask0)
    pos = 0
    mask = ""
    blocks0 = []
    for pp in p:
        blocks0.append(s[pp.start():pp.end()])
        preblock = s[pos:pp.start()]
        mask += preblock + "!" * (pp.end()-pp.start())
        pos = pp.end()
    mask += mask0[pos:]
    if pos == 0:
      break
    else: postblock = s[pos:]
    blocks = blocks0
    mask0 = mask
  block = None
  if len(blocks) == 1: 
    block = blocks[0][1:-1]
  if len(blocks) > 1:
    raise Exception("compile error: invalid statement\n%s\nStatement must contain a single {} block" % s)
    if len(postblock.strip()) != 0:
      raise Exception("compile error: invalid statement\n%s\nStatement must be empty after {} block" % s)
  preblock = preblock.strip()
  
  blockcomment = ""
  
  p = quotematch.finditer(preblock)
  pos = 0
  mask = ""
  for pp in p:
      mask += preblock[pos:pp.start()] + "*" * (pp.end()-pp.start())
      pos = pp.end()
  mask += preblock[pos:]  
  comment = mask.find("#")
  if comment > -1:     
    blockcomment = preblock[comment+1:].strip('\n') + '\n'
    preblock = preblock[:comment]
  
  if block != None:
    currblock = block
    while 1:
      len0 = len(currblock)
      currblock = currblock.lstrip().lstrip("\n")
      if len(currblock) == len0: break
      
    pp = triplequotematch.search(currblock)
    if pp != None and pp.start() == 0:
      blockcomment += currblock[pp.start()+len('"""'):pp.end()-len('"""')]
            
  name = None
  title = None
  if len(preblock) > 0:
    name = preblock.split()[0]
    title = preblock[len(name):].strip()
    
  return name,title,block, blockcomment

def includetype(classname, inherited, textblocks, typedefs):
  if validvar2(inherited) == False:
    raise Exception("compile error: invalid type definition\n%s\nCannot inherit from non-Spyder class '%s'" % (l, inherited))
  if inherited not in typedefs:
    import Spyder
    getattr(Spyder, inherited)      
  if inherited not in typedefs:
    raise Exception("compile error: Spyder class %s defines from undefined Spyder class %s" % (classname, inherited))    
  parenttype = typedefs[inherited]    
  for v in parenttype:
    if not v.startswith("constructor_"):
      if v not in textblocks:
        textblocks[v] = copy(parenttype[v])
      elif v == "members":
        textblocks[v] += [tuple(list(i) + ["include"]) for i in parenttype[v]]
      else:
        textblocks[v] += parenttype[v]


def typedef(l, title, block, typedefs, currmodule, currfile):
  from .define import unimethodlist, unimethods, macros
  from . import silent
  
  neededtypes = set()
  if len(title) == 0 or block == None:
    raise Exception("compile error: invalid type definition\n%s\nType definition must have a proper title" % l)
  t = [tt.strip() for tt in title.split(":")]
  if len(t) > 2:
    raise Exception("compile error: invalid type definition\n%s\nInvalid title: %s" % (l, title))
  typename = t[0]
  if  validvar(typename) == False:
    raise Exception("compile error: invalid type definition\n%s\nInvalid title: %s" % (l, title))  
  if typename in reservedtypes:
    raise Exception("compile error: reserved type name %s" % (typename))
  if typename in typedefs:
    raise Exception("compile error: type name %s has already been defined" % (typename))
  oldtextblocks = {}
  oldtextblocks["members"] = []
  oldtextblocks["deleted_members"] = []
  inheriteds = []
  if len(t) == 2:
    inheriteds = [inherited.strip() for inherited in t[1].split(",")]
    for inherited in inheriteds:
      includetype(title, inherited, oldtextblocks,typedefs)      
  
  textblocks = oldtextblocks
  inherited_array = False
  for inherited in inheriteds:
    inh2 = inherited    
    if inherited.endswith("Array"): 
      inherited_array = True
      while inh2.endswith("Array"): inh2 = inh2[:-len("Array")]
      assert len(inheriteds) == 1 #only single inheritance from Array
    neededtypes.add(inh2)  
  if inherited_array:
    if "validate" not in textblocks: textblocks["validate"] = ""
    textblocks["validate"] += "\nself.__arrayvalidate__()\n"
        
  filtersource = "Spyder-generated class\n\n module %s\n file %s\n\n" % (currmodule, repr(currfile))
  
  filtersource0, commentblock = typedefparse(block, textblocks, macros, typedefs, neededtypes) 

  #textblocks gets updated...
  
  filtersource += commentblock + '\n'
  
  filtersource += 70 * "-" + "\nSpyder definition:\n\n"          
  if len(inheriteds):
    filtersource += "Type %s(%s) {\n" % (typename, ",".join(inheriteds))
  else: 
    filtersource += "Type %s {\n" % (typename)  

  filtersource += filtersource0
  textblocks["members"] = [a for a in textblocks["members"] if a[0] != "Delete"]
  filtersource += "}\n"
  filtersource += '\n' + 70 * "-" + '\n'
  
  code = ""  
  unimethodlist0 = copy(unimethodlist)
  unimethodlist0.remove("__class")
  unimethodlist0.remove("__endclass")
  endclass = {}
  allmethods = list(["__class",] + unimethodlist0)
  members = textblocks["members"]
  deleted_members = textblocks["deleted_members"]
  for u in allmethods:
    block = None
    if u in textblocks: block = textblocks[u]
    ok = True
    try:
      unimethodcode, endclasscode = unimethods[u](typename, inheriteds, filtersource, members, deleted_members, block)
    except:
      import traceback  
      print("*" * 60)
      traceback.print_exc()
      print("*" * 60)
      ok = False
    if not ok: raise ImportError
    endclass[u] = endclasscode
    if unimethodcode != None: code += unimethodcode
  block = None
  if "__endclass" in textblocks: block = textblocks["__endclass"]    
  endclasscode = unimethods["__endclass"](typename, inheriteds, filtersource, members, deleted_members, block, endclass)
  if endclasscode != None: code += endclasscode
  #safe_eval(code, True) DISABLED! unimethods must check!
  if "__endclass" in textblocks:
     testcode = "class bogus:\n" + textblocks["__endclass"]
     safe_eval(testcode)
  typedefs[typename] = textblocks
  code += "\n"
  return code, neededtypes

def typedefparse(block, textblocks, macros, typedefs, neededtypes):  
  filtersource = ""
  commentblock = ""
    
  lines = getlines(block)
  defmode = False
  currindent = 0
  prelines = 0
  memberblock = None
  while len(lines) > 0:
    if prelines > 0: prelines -= 1
    l = lines[0].strip()
    l2 = lines[0].replace('\t', "  ")
    lines = lines[1:]
    if len(l) == 0: continue
    if l.startswith("##"):
      l = l[2:].lstrip()
      ll = l.split()
      name = ll[0]
      if name not in textblocks: textblocks[name] = ""
      textblocks[name] += l[len(name)+1:] + '\n'
      continue
    if l.find("#") == 0: 
      spacing = l2[:len(l2)-len(l)]
      filtersource += spacing + l + '\n'
      continue
    if defmode == False and l.startswith('"""'):
      spacing = l2[:len(l2)-len(l)]
      l = l[3:]
      l = '\n' + l[:l.index('"""')]
      filtersource += spacing + "#!#!#!\n" + spacing + l.strip('\n').strip() + '\n' + spacing + "#!#!#!\n"
      if memberblock: textblocks["members"][-1][2] += l[len('\n'):]
      elif memberblock == None: commentblock += l[len('\n'):]
      continue
    if defmode == True:
      indent = len(l2) - len(l2.lstrip())
      if indent == currindent:
        defmode = False
    if defmode == False:
      if l2.lstrip().startswith("def "):
        currindent = len(l2) - len(l2.lstrip())
        defmode = True
        if prelines == 0: filtersource += l2 + '\n'
    if defmode == True or l2.lstrip().startswith("@"):      
      if "__endclass" not in textblocks: textblocks["__endclass"] = ""
      textblocks["__endclass"] += l2 + '\n'
      continue
    if prelines == 0: filtersource += l2 + '\n'
    
    name,title,block,blockcomment = getblock(l)
    memberblock = False
    if block != None:
      if title != "" and title != None:
        raise Exception("compile error: unimethod statement must be <name> {...}\n%s" % (l))
      if name not in textblocks: textblocks[name] = ""
      bb = block.split('\n')
      spaces = -1
      for l in bb:
        if len(l.strip()) == 0: continue
        if spaces == -1:
          spaces = len(l) - len(l.lstrip())
        textblocks[name] += l.rstrip('\n')[spaces:] + '\n'
    else:
      triggermacro = False
      for m in macros:
        mm = m(name, title)
        if mm != None:
          triggermacro = True
          oldlen = len(lines)
          lines = getlines(mm) + lines
          prelines += len(lines) - oldlen
          break
      if triggermacro == False:
        if name == "Delete":
          if title in textblocks:
            textblocks.pop(title)  
          else:
            textblocks["deleted_members"].append(title)
          textblocks["members"] = [a for a in textblocks["members"] if a[1].split("=")[0].rstrip() != title]
          textblocks["members"].append([name,title, blockcomment])
        elif name == "Include":
          includetype(name, title, textblocks,typedefs)
        elif not validvar2(name):
          raise TypeError("Invalid member name %s" % name)
        else:
          neededtypes.add(name)
          memberblock = True 
          textblocks["members"].append([name,title, blockcomment])
          
  members0 = {}
  memberlist = set()
  counter = 0
  for tb in textblocks["members"]:
    membertype, title, blockcomment = tb[:3]
    if membertype == "Delete": 
      members0[len(members0)] = [membertype, title, blockcomment, len(members0)]
      continue
    membername = title.split("=")[0].rstrip()
    if membername.replace("_","x").isalnum() == False:
      raise Exception("compile error: invalid membername %s" % (membername))
    if membername in reserved_membernames:
      raise Exception("compile error: reserved membername %s" % (membername))    
    if membername in members0 and members0[membername][0] != None:
      if len(tb) > 3 and tb[3] == "include": continue
      raise Exception("compile error: duplicate membername %s" % membername)
    if membername in members0: members0[membername][0] = membertype
    else: members0[membername] = [membertype, title, blockcomment, len(members0)]
  members00 = {}
  for m in members0:
    if members0[m][0] != None:
      members00[members0[m][3]] = (members0[m][0], members0[m][1], members0[m][2])
  memnrs = copy(list(members00.keys()))
  memnrs.sort()
  textblocks["members"] = []
  for nr in memnrs: textblocks["members"].append(members00[nr])
  #print "MEMBERS", textblocks["members"]
  return filtersource, commentblock

def statement(l0, typedefs,currmodule,currfile):
  l = l0.strip()
  if len(l) == 0: return "", set(), set()
  if l[0] == "#": return "", set(), set()

  name,title,block,blockcomment = getblock(l)
  
  if name == "Type":
    code,neededtypes = typedef(l,title, block, typedefs,currmodule,currfile)
    return code, set(), neededtypes
  from .define import statements
  for statement in statements:
    if name == statement:
      supportedtypes = set()
      r, sup = statements[statement](l, title, block)
      for suptyp in sup:
        while suptyp.endswith("Array"):
          suptyp=suptyp[:-len("Array")]
        if suptyp.startswith("Resource"): suptyp=suptyp[len("Resource"):]  
        if suptyp != "None": supportedtypes.add(suptyp)
      return r, supportedtypes, set()
  return l0 + '\n', set(), set()

def compile_obj(code, name, sourcefile, bytecode, recompile):
  bytecode_ok = False  
  if recompile is not True and bytecode is not None:
    try:
      magic = bytecode[:4]
      assert magic == imp.get_magic()  
      if recompile is not False:
        timestamp = struct.unpack("<l",bytecode[4:8])[0]
        lastmodified = os.stat(sourcefile)[9]
        assert (timestamp >= lastmodified), (timestamp, lastmodified)
      header = 8  
      pyver = sys.version_info[:2]
      if pyver[0] == 3 and pyver[1] >= 3: #python >= 3.3
        header = 12
      obj = marshal.loads(bytecode[header:])
      bytecode_ok = True
    except:
      pass  
  if not bytecode_ok:
    obj = compile(code, name, "exec")
    try:
      import importlib
      importlib.invalidate_caches()
    except:
      pass
  return obj, bytecode_ok

def visit_node(node, defined, used):
  nodetype = node.__class__.__name__
  if nodetype == "Name":
    name = list(ast.iter_fields(node))[0][1]
    if validvar2(name):      
      used.add(name)
  elif nodetype == "ClassDef":
    name = list(ast.iter_fields(node))[0][1]
    if validvar2(name):      
      defined.add(name)
  children = ast.iter_child_nodes(node)
  for child in children:
    visit_node(child, defined, used)

def update_dependencies(definedtypes, neededtypes, usedtypes):
  def chop(t):
    tt = t
    while tt.endswith("Array"): tt = tt[:-len("Array")]
    if tt.startswith("Resource"): tt = tt[len("Resource"):]    
    return tt
  for t in list(neededtypes):
    if t in definedtypes: neededtypes.remove(t)
    tt = chop(t)
    if t != tt:
      if t in neededtypes: neededtypes.remove(t)
      if tt not in neededtypes and tt not in definedtypes: neededtypes.add(tt)
  for t in list(usedtypes):
    tt = chop(t)
    if tt in neededtypes or tt in definedtypes: usedtypes.remove(t)

def get_called_types(filename, code):
  astr = ast.parse(code, filename)
  defined, used = set(), set()
  visit_node(astr, defined, used)
  return defined, used
  
def write_types(filename, definedtypes, neededtypes, usedtypes):
  try:
    f = open(filename, "w")
  except IOError:
    return
  for t in definedtypes:
    print("DEFINED",t,file=f)
  for t in neededtypes:
    print("NEEDED",t,file=f)
  for t in usedtypes:
    print("USED",t,file=f)
  f.close()

from .loader import file_load, file_exists
def get_types(filename, code, typef, fresh, cache=False):
  ok = False
  if file_exists(typef) and fresh:
    try:
      defined, needed, used = set(),set(),set()
      f = file_load(typef)
      types = f.readlines()
      f.close()

      for l in types:
        mode,typ = l.split()
        if mode == "DEFINED": defined.add(typ)
        elif mode == "NEEDED": needed.add(typ)
        elif mode == "USED": used.add(typ)
      ok = True
    except:
      ok = False
  if not ok:
    needed = None
    defined, used = get_called_types(filename, code)
    if cache:
      write_types(typef, defined, set(), used)    
  return defined, needed, used
 

def spydercompile(s,tempdir,recompile,temp_prefix,corechecksum,typedefs,currmodule,currfile):
  from . import __types__, spyderdir

  tempdir = tempdir.replace("#SPYDERDIR", spyderdir)
  tempdir = tempdir.replace("\\", os.sep).replace("/", os.sep)
  
  cachename = os.path.split(currfile)[1]
  
  definedtypes = set()
  neededtypes = set()
  usedtypes = set()
  
  cachepat = tempdir + os.sep + temp_prefix + "." + cachename
  cachefile = cachepat + ".py"
  cachefile2 = cache_from_source(cachefile)
  cachefile_types = cachepat + ".TYPES"
  cachefile_typedefs = cachepat + ".TYPEDEFS"

  codesum = None  
  if file_exists(cachefile) and file_exists(cachefile_typedefs):
    cacheftd = file_load(cachefile_typedefs)
    cache_typedefs = [l.strip() for l in cacheftd.readlines()]
    cacheftd.close()
    
    cachecoresum = int(cache_typedefs[0])
    cachecodesum = int(cache_typedefs[1])

    fresh = True 
    if recompile == True:
      fresh = False
    elif recompile != False and cachecoresum != corechecksum:
      fresh = False
    else:
      cachef = file_load(cachefile)
      cachecode = cachef.read()
      cachef.close()

      if recompile is None:
        codesum = zlib.adler32(s.encode('UTF-8')) & 0xffffffff
        if codesum != cachecodesum: fresh = False
         
    if fresh:
      definedtypes, neededtypes, usedtypes = get_types(cachefile, cachecode, cachefile_types, True)

    if fresh:  
      lnr = 2
      while lnr < len(cache_typedefs) - 1:
        t = cache_typedefs[lnr]
        deff = cache_typedefs[lnr+1]
        typedefs[t] = eval(deff)
        lnr += 2
      #print("No recompile %s" %  currmodule+"::"+currfile, len(cachecode))

      bytecode_ok = False
      cachecode2 = None
      if file_exists(cachefile2):
        cachef2 = file_load(cachefile2,"rb")
        cachecode2 = cachef2.read()
        cachef2.close()
      obj, bytecode_ok = compile_obj(cachecode,currmodule+"::"+currfile,currfile,cachecode2,recompile)
      if not bytecode_ok:
        try:
          os.remove(cachefile2)
        except:
          pass
        import py_compile
        try:
          py_compile.compile(cachefile,cachefile2,currmodule+"::"+currfile)
        except:
          pass          
      for t in definedtypes:
        if t not in typedefs: typedefs[t] = []
      return cachecode, obj, definedtypes, neededtypes, usedtypes

  from . import silent
  if not silent:
    print("Recompile %s" %  currmodule+"::"+currfile)
  if codesum is None: codesum = zlib.adler32(s.encode('UTF-8')) & 0xffffffff
  
  code = ""
  
  lines = getlines(s)
  for l in lines:
    code0,supportedtypes,needed = statement(l, typedefs, currmodule, currfile)
    usedtypes.update(supportedtypes)
    neededtypes.update(needed)
    if code0 is not None: code += code0
    
  try:
    definedtypes, dmmy, used = get_types(cachefile, code, cachefile_types, fresh=False)
  except SyntaxError:
    cachef = open(cachefile, "w")
    cachef.write(code)
    cachef.close()
    raise  
  usedtypes.update(used)
  for t in definedtypes:
    if t not in typedefs: typedefs[t] = []
  
  cachef = open(cachefile, "w")
  cachef.write(code)
  cachef.close()
  try:
    import importlib
    importlib.invalidate_caches()
  except:
    pass  
  
  import py_compile
  try:
    py_compile.compile(cachefile,cachefile2,currmodule+"::"+currfile)
  except:
    pass
  update_dependencies(definedtypes, neededtypes, usedtypes)
  write_types(cachefile_types, definedtypes, neededtypes, usedtypes)
  cacheftd = open(cachefile_typedefs, "w")
  print(corechecksum,file=cacheftd)
  print(codesum,file=cacheftd)  
  for t in definedtypes:
    print(t,file=cacheftd)
    print(repr(typedefs[t]),file=cacheftd)
  cacheftd.close()

  bytecode_ok = False
  cachecode2 = None
  if file_exists(cachefile2):
    cachef2 = file_load(cachefile2,"rb")
    cachecode2 = cachef2.read()
    cachef2.close()    
  obj, bytecode_ok = compile_obj(code,currmodule+"::"+currfile,currfile,cachecode2,recompile)
  if not bytecode_ok:
    try:
      os.remove(cachefile2)
    except:
      pass    
  return code, obj, definedtypes, neededtypes, usedtypes
 
