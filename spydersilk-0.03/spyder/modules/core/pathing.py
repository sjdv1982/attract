# Copyright 2010-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

# Assumptions of the conversion engine:
# - Last defined, first tried
# - Converters do not have side effects, they do not change outside variables (but methods do!)
# - Converters do not change the variable they operate on
# - A converter can fail by returning None, indicating that the converter cannot operate on the provided values; this allows specialization
# - A converter can be dynamically disabled and enabled
# - The path to a type doesn't matter, 
#    i.e. the conversion A to B has a unambiguous result independent of the conversion path
#    This also means that every Spyder type is a point in a path that is visited only once
# - There are automatic SPLIT converters from XArray to X, and from X to an XArray consisting of a single X
#    However, SPLIT converters only work for method search
#     and the conversion engine will NOT automatically convert XArray to YArray if there is a converter X=>Y

# Split None converters
# - They can only be used with method finding paths, not conversion paths
# - They are considered the endpoint of a partial path: 
#    After the split, splitnonecounter is incremented and the method finding is started all over for the result
#   (including a reset of visited and failed)
# - A conversion path can have at maximum "max_splitnone" SPLIT-None converters

from functools import partial

def validvar(s):
  if s.replace("_", "x").isalnum() == False: return False
  if s[0].isupper() == False: return False
  if len(s) > 1 and s == s.upper(): return False
  return True

class PathList:
  def __init__(self):
    self.converterindices  = list(range(1,len(converters)+1))
    self.possible = set()
    self.busy = False
      
pathlists = {}
methodpathlists = {}
max_splitnone = 20

class spyderconverter:
  def __init__(self, intype, outtype, pointer):
    self.intype = intype
    self.outtype = outtype
    self.pointer = pointer
    self.enabled = True
  def enable(self):
    self.enabled = True
  def disable(self):
    self.enabled = False

class spydermethod:
  def __init__(self, pointer):
    self.pointer = pointer
    self.enabled = True
  def enable(self):
    self.enabled = True
  def disable(self):
    self.enabled = False

converters = []
methods = {} #dictionary that records the methods of each type
methodtypes = {} #dictionary that records the types of each method

define_functioncounter = 1
conversioncounter = 0

def defineconverter(intype, outtype, pointer):
  global pathlists, methodpathlists
  pathlists, methodpathlists = {}, {}
  if not validvar(intype):
    raise Exception("compile error: type %s is not a valid Spyder type" % intype)
  if not validvar(outtype):
    if outtype != "None" or pointer != "SPLIT":
      raise Exception("compile error: type %s is not a valid Spyder type" % outtype)
  c = spyderconverter(intype, outtype, pointer)
  converters.append(c)
  return c

def definemethod(method,type,pointer):
  global pathlists, methodpathlists, methods, methodtypes
  pathlists, methodpathlists = {}, {}
  if not validvar(type):
    raise Exception("compile error: type %s is not a valid Spyder type" % type)
  if type not in methods: methods[type] = {}
  """
  if method in methods[type]:
    raise Exception("compile error: method \"%s\" for type %s has already been defined" % (method,type))
  #print "METHOD DEFINED", method, type
  methods[type][method] = pointer
  """
  if method not in methods[type]: methods[type][method] = []
  newmethod = spydermethod(pointer)
  methods[type][method].append(newmethod)

  if method not in methodtypes: methodtypes[method] = set()
  methodtypes[method].add(type)
  return newmethod

def callall(pointers, *a, **aa):
  ret = []
  for p in pointers:
    ret.append(p(*a, **aa))
  return ret

def execute_path(i,arg,path,failed,visited,method,splitcounter,methodindex):
  #This will execute the path that was found previously by find_convert or find_method
  assert i == arg.typename()  
     
  if len(path) == 0: #We have reached the end of the path!
    arg.__conversionstack__ = list(__conversionstack__.l)
    if method == None:
      #In case of conversions, return the Spyder object of the requested type      
      return arg      
    else:
      #In case of method finding, return a callable object 
      return partial(methods[i][method][methodindex].pointer, arg)
      
  cnr = path[0]
  c = converters[cnr-1]
  f = c.pointer
  pathtail = path[1:]

  if f == "SPLIT":     
    assert method != None #SPLIT converters can only be used with method finding
    if c.outtype == "None":  
      assert len(path) == 1 #our path is partial; it must be the last converter
      if splitcounter == max_splitnone: return None
      splitcounter += 1                  
      result = []
      for anr,a in enumerate(arg):
        __conversionstack__.l.append(anr)
        #print "SPLIT NONE", anr+1, __conversionstack__.l
        stackbackup = list(__conversionstack__.l)        
        r = do_method(a.typename(),method,a,splitcounter)
        __conversionstack__.l = stackbackup
        if r == None:
          result = None
          break        
        else: result.append(r)
      if result != None: return partial(callall,result)
    else:      
      result = []
      for a in arg:
        stackbackup = list(__conversionstack__.l)
        r = execute_path(c.outtype, a, pathtail, failed,visited,method,splitcounter,methodindex)
        __conversionstack__.l = stackbackup
        if r == None:
          result = None
          break        
        else: result.append(r)
      if result != None: return partial(callall,result)          
  elif f == "CAST":
    result = spyder.__types__[c.outtype](arg)
  else: 
    assert c.intype == arg.typename()
    result = f(arg)

  if result == None:     
    path_ok = False
    failed.add(cnr)
  else:  
    path_ok = True
    if not hasattr(result, "typename") or result.typename() != c.outtype:
      result = spyder.__types__[c.outtype](result)    
    
  while 1:  
    if path_ok:
      visited.add(i)
      stackbackup = list(__conversionstack__.l)              
      arg = result      
      result = execute_path(c.outtype, arg, pathtail, failed, visited, method,splitcounter,methodindex)
      __conversionstack__.l = stackbackup      
    else:
      if method != None:
        newpath,methodindex = find_method(i, method, failed, visited)      
      else:
        o = converters[path[-1]-1].outtype      
        newpath = find_convert(i,o,failed,visited)
        methodindex = None
      if newpath == None: 
        return None
      stackbackup = list(__conversionstack__.l)        
      result = execute_path(i,arg,newpath,failed,visited,method,splitcounter,methodindex)
      __conversionstack__.l = stackbackup
      
    if result != None: return result    
    path_ok = False
    
    
  

def find_convert(i,o,failed,visited):
  #tries to find a conversion path from Spyder type i to Spyder type o
  # given a set of failed converter indices
  #  and a visited list of Spyder types

  #if (i,o) not in pathlists:
  #  pathlists[(i,o)] = PathList()
  #p = pathlists[(i,o)]
  p = pathlists.setdefault((i,o), PathList())

  pos = len(p.converterindices)
  newvisited = set(visited); newvisited.add(i)
  while pos > 0:
    pos -= 1
    cnr = p.converterindices[pos]    
    c = converters[cnr-1]    
    
    #Prune:
    #- Converters that do not operate on the current type
    #- Converters that lead nowhere, even without visited types or failed/disabled converters
    prune = True
    if cnr in p.possible: prune = False
    elif c.intype == i and c.pointer != "SPLIT" and c.outtype not in visited:
      old_busy = p.busy
      p.busy = True
      if c.outtype == o or find_convert(c.outtype,o,set(),newvisited) != None: 
        p.possible.add(cnr)
        prune = False
      p.busy = old_busy
    if prune:
      p.converterindices.pop(pos)
      continue    
    
    #Ignore:
    #- Disabled converters
    #- Failed converters
    #- Converters that lead to visited types
    if not c.enabled: continue
    if cnr in failed: continue
    if c.outtype in visited: continue
    
    #Let's see where this one leads...
    if c.outtype == o: return [cnr]
    else: 
      #newvisited = set(visited)
      #newvisited.add(i)
      ret = find_convert(c.outtype,o,failed,newvisited)
      if ret != None: return [cnr] + ret
      
      #It led nowhere, just continue...
      
  return None
  
def do_convert(i, o, arg):
  failed = set()
  visited = set()
  path = find_convert(i,o,failed,visited)
  if path == None:
    raise Exception("No conversion path from %s to %s" % (i,o))
  conv = execute_path(i,arg,path,failed,visited,None,0,None)
  if conv == None: 
    raise Exception("Failed to convert %s to %s" % (i,o))
  return conv
    
def convert(intype, outtype, arg,deepcopy):
  stackbackup = list(__conversionstack__.l)
  i = intype.typename()
  o = outtype.typename()
  if i == o:     
    arg.__conversionstack__ = list(__conversionstack__.l)
    ret = arg
  else:
    ret = do_convert(i,o,arg)  
  __conversionstack__.l = stackbackup
  if not deepcopy: return ret
  return type(ret)(ret)


#####


def find_method(i,m,failed,visited):
  #tries to acquire through conversion a method m for Spyder type i
  # given a set of failed converter indices
  #  and a visited list of Spyder types
  
  #if (i,m) not in methodpathlists:
  #  methodpathlists[(i,m)] = PathList()
  #p = methodpathlists[(i,m)]
  p = methodpathlists.setdefault((i,m), PathList())
  if p.busy: return

  mt = methodtypes[m]
  pos = len(p.converterindices)
  newvisited = set(visited); newvisited.add(i)
  while pos > 0:
    pos -= 1
    cnr = p.converterindices[pos]    
    c = converters[cnr-1]    
    
    #Prune:
    #- Converters that do not operate on the current type
    #- Converters that lead nowhere, even without failed/disabled converters
    prune = True
    if cnr in p.possible: prune = False    
    elif c.intype == i:
      if c.outtype not in visited:
        old_busy = p.busy
        p.busy = True      
        if c.outtype in mt or c.outtype == "None" or find_method(c.outtype,m,set(),newvisited) != None: 
          p.possible.add(cnr)
          prune = False
        p.busy = old_busy
    if prune:
      p.converterindices.pop(pos)
      continue    

    #Ignore:
    #- Disabled converters
    #- Failed converters
    #- Converters that lead to visited types
    if not c.enabled: continue
    if cnr in failed: continue
    if c.outtype in visited: continue
    
    #Let's see where this one leads...
    if c.outtype in mt:
      #if even a single method is enabled, it is OK
      for mnr,meth in enumerate(methods[c.outtype][m]):
        if meth.enabled: return [cnr], mnr
    elif c.outtype=="None": #the latter is a SPLIT-None converter
      return [cnr], None
    else: 
      #newvisited = set(visited)
      #newvisited.add(i)
      ret = find_method(c.outtype,m,failed,newvisited)
      if ret != None: 
        return [cnr] + ret[0], ret[1]
      
      #It led nowhere, just continue...
      
  return None
  
def do_method(i,m,arg,splitcounter):
  global methods
  if i not in methods: methods[i] = {} 
  if m not in methodtypes: methodtypes[m] = {} 

  if m in methods[i]:
    for mm in methods[i][m]:
      if not mm.enabled: continue
      arg.__conversionstack__ = list(__conversionstack__.l)
      ret = partial(mm.pointer, arg)
      return ret

  failed = set()
  visited = set()
  path = None
  if len(methodtypes[m]):
    path = find_method(i,m,failed,visited)
  if path == None:
    raise AttributeError("No path to acquire method %s for %s" % (m,i))
  path,methodindex = path
  ret = execute_path(i,arg,path,failed,visited,m,splitcounter,methodindex)
  if ret == None:
    raise AttributeError("Failed to acquire method %s for %s" % (m,i))
  return ret
  
def method(type, methodname, arg):  
  stackbackup = list(__conversionstack__.l)
  global conversioncounter
  conversioncounter += 1
  try:
    t = type.typename()
    ret = do_method(t,methodname,arg,0)
  finally:  
    conversioncounter -= 1
    if conversioncounter == 0: 
      spyder.core.__clear_temp()
    __conversionstack__.l = stackbackup
  return ret
  
#the temporary dictionary is reset after every invocation of the conversion engine  
  
__spyder_temp = {}
def __set_temp(key, value):
  global __spyder_temp
  __spyder_temp[key] = value
  
def __get_temp(key):
  global __spyder_temp
  return __spyder_temp.get(key,None)

def __clear_temp():
  global __spyder_temp
  __spyder_temp = {}

  
