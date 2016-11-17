from . import namestr, log

from collections import OrderedDict
import sys, os, weakref
import functools, operator
import numpy as np
import traceback

class TVDepsTask(object):
  def __init__(self, task):
    self.tvtask = weakref.ref(task)
    self.data = []
    self.resources = []
    self._cachekeys = []
  def _cachekeys_load(self):
    for a in self._cachekeys:
      arr = np.load(a._cache.cachekey)    
      a.set_data(arr)
      a._cache.cachekey = None
  def cleanup(self, graph):
    self.data = [d for d in self.data if d in graph.data]
    self.resources = [d for d in self.resources if d in graph.resources]
  def empty(self):
    return (not self.data) and (not self.resources)
  def __str__(self):
    s = "    " + str(self.tvtask()) + "\n"    
    namestrs = [namestr(d) for d in self.data]
    s += "    Data dependencies:\n      " + " ".join(namestrs) + "\n"
    if len(self.resources):
      namestrs = [namestr(r) for r in self.resources]
      s += "    Resource dependencies:\n      " + " ".join(namestrs) + "\n"
    return s
  def dict(self, state):
    t = self.tvtask()
    tt = state.find_task(t)
    return {"tvtask": tt, "data":self.data,"resources":self.resources}
  @classmethod
  def fromdict(cls, dic, state):
    self = object.__new__(cls)
    self.tvtask = weakref.ref(state.tasks[dic["tvtask"]])
    self.data = dic["data"]
    self.resources = dic["resources"]
    self._cachekeys = []
    return self
  
class TVDepsResource(object):
  def __init__(self, tvdata):
    self.tvdata = lambda: None
    if tvdata is not None:
      if isinstance(tvdata, weakref.ref):
        self.tvdata = tvdata
      else:
        self.tvdata = weakref.ref(tvdata)
    self.data = []
    self.resources = []
    self.allocators = []
  def cleanup(self, graph):
    self.data = [d for d in self.data if d in graph.data]
    self.resources = [d for d in self.resources if d in graph.resources]
    self.allocators = [d for d in self.allocators if d in graph.allocators]
  def empty(self):
    return (not self.data) and (not self.resources) and (not self.allocators)
  def __str__(self):
    s = ""    
    if len(self.data):
      namestrs = [namestr(d) for d in self.data]
      s += "    Data dependencies:\n      " + " ".join(namestrs) + "\n"
    if len(self.resources):
      namestrs = [namestr(d) for d in self.resources]
      s += "    Resource dependencies:\n      " + " ".join(namestrs) + "\n"
    if len(self.allocators):
      namestrs = [namestr(d) for d in self.allocators]
      s += "    Allocators:\n      " + " ".join(namestrs) + "\n"
    return s
  def dict(self, state):
    dic = {"data":self.data,"resources":self.resources, "allocators": self.allocators}
    t = self.tvdata()
    if t is not None:
      if isinstance(t, TVArray):
        dic["tvdatatype"] = "array"
        tt = state.find_data(t)
      else:
        dic["tvdatatype"] = "container"
        tt = state.find_data_container(t)        
      dic["tvdata"] = tt
    return dic 
  @classmethod
  def fromdict(cls, dic, state):
    self = object.__new__(cls)
    if "tvdata" in dic:
      if dic["tvdatatype"] == "array": 
        d = state.data[dic["tvdata"]]
      elif dic["tvdatatype"] == "container":
        d = state.data_containers[dic["tvdata"]]
      else:
        raise TypeError(dic["tvdatatype"])
      self.tvdata = weakref.ref(d)
    else:
      self.tvdata = lambda: None
    self.data = dic["data"]
    self.resources = dic["resources"]
    self.allocators = dic["allocators"]
    return self

def make_hashtree(dataname, graph):
  d = graph.data[dataname]
  data = d.tvdata()
  assert data is not None and data._cache_active == False, dataname
  t = d.task
  assert t is not None, dataname
  task = graph.tasks[t].tvtask()
  subtrees = []
  for inr, i in enumerate(task._inputs):
    if isinstance(i, TVArray):       
      ddname = i.name()
      #print "HASHING", namestr(ddname), task
      ddname2 = ddname
      if isinstance(i, TVReadShard) and not i._parent._cache_active:
        ddname2 = ddname[:2]      
      dd = graph.data[ddname2]
      if not i._cache_active and (not hasattr(i, "_hash") or i._hash is None):        
        if dd.hashtree is None: 
          dd.hashtree = make_hashtree(ddname2, graph)
        subtree = ("hashtree", dd.hashtree)
      else:
        subtree = ("hash", ddname)
      if isinstance(i, TVReadShard) and not i._parent._cache_active:
        subtree = ("shard", (i._offset, i.shape, subtree))
    else:
      subtree = ("signature", task.signature[inr])
    
    subtrees.append(subtree)
    
  return HashTree((task.action.name(), tuple(subtrees)))

class HashTree(tuple):
  bound = None
  def bindable(self, graph):
    taskname, subtrees = self
    for subtree in subtrees:
      mode, sub = subtree
      shard = False
      if mode == "shard":
        shard = True
        offset, shape, subtree = sub
        mode, sub = subtree
      if mode == "hash": #standard hash
        d = graph.data[sub]
        if d.task is not None:
          i = d.tvdata()
          if i.get_hash() is None:
            return False
      elif mode == "hashtree":
        assert isinstance(sub, HashTree)
        ok = sub.bindable(graph)
        if not ok: return False
      elif mode == "signature":
        pass
    return True  
  def bind(self, graph):
    if self.bound is not None: 
      return self.bound
    taskname, subtrees = self
    ret = [taskname]
    for subtree in subtrees:
      mode, sub = subtree
      shard = False
      if mode == "shard":
        shard = True
        offset, shape, subtree = sub
        mode, sub = subtree
      if mode == "hash": #standard hash
        d = graph.data_unpruned[sub] #unpruned
        i = d.tvdata()
        hash = i.get_hash() 
        assert hash is not None, sub
        append = hash        
      elif mode == "hashtree":
        assert isinstance(sub, HashTree)
        append = sub.bind(graph)
      elif mode == "signature":
        append = sub
      
      if shard:
        append = (offset, shape, append)
      ret.append(append)    
    
    ret = tuple(ret)    
    self.bound = ret
    return ret
    
    
class TVDepsData(object):
  hashtree = None
  def __init__(self, tvdata):
    self.tvdata = lambda: None
    if tvdata is not None:
      self.tvdata = weakref.ref(tvdata)
    self.task = None
    self.resource = None
  def cleanup(self, graph):
    if self.task is not None and self.task not in graph.tasks:
      self.task = None
    if self.resource is not None and self.resource not in graph.resources:
      self.resource = None      
  def empty(self):
    return self.task is None and self.resource is None    
  def __str__(self):
    if self.task is None and self.resource is None: return ""
    if self.resource is not None:
      s = "    Resource dependency:\n      " + str(namestr(self.resource)) + "\n"
    if self.task is not None:
      s = "    Task dependency:\n      " + str(self.task) + "\n"
    return s
  def dict(self, state):
    dic = {}
    t = self.tvdata()
    if t is not None:
      tt = state.find_data(t)
      dic["tvdata"] = tt
    if self.task is not None:
      dic["task"] = self.task
    if self.resource is not None:
      dic["resource"] = self.resource
    if self.hashtree is not None:
      self.hashtree.bound = None
      dic["hashtree"] = self.hashtree
    return dic 
  @classmethod
  def fromdict(cls, dic, state):
    self = object.__new__(cls)
    if "tvdata" in dic:
      self.tvdata = weakref.ref(state.data[dic["tvdata"]])
    else:
      self.tvdata = lambda: None
    self.task = dic.get("task", None)
    self.resource = dic.get("resource", None)
    self.hashtree = dic.get("hashtree", None)
    return self

class TVDepsGraph(object):
  cachedir = None
  cache_indices = None
  def __init__(self, state, clone=False):
    self.tasks = OrderedDict()  
    self.data = OrderedDict()
    self.resources = OrderedDict()
    self.allocators = OrderedDict()
    self._evaluated = False   
    self.state = state
    
    if clone: return
    
    from .TVAction import TVAction, TVTask
    from .TVArray import TVArrayBase, TVArray, TVBaseShard, TVReadShard, TVArraySeries  
  
    #1. Initial build of depsgraph from parent register
    for a in self.state.tasks:
      task = TVDepsTask(a)
      self.tasks[len(self.tasks)+1] = task
    for a in self.state.data:  
      assert isinstance(a, TVArray)
      
      aname = a.name()
      d = TVDepsData(a)
      assert aname not in self.data, namestr(aname)
      self.data[aname] = d 
      a._cache.graph = self
      if isinstance(a, TVBaseShard): continue
    
      aname2 = aname
      if a._parent is not None and aname[-1] is None:
        aname2 = aname[0], 0 
      assert aname2 not in self.resources, namestr(aname2)      
      r = TVDepsResource(a)
      self.resources[aname2] = r
    for a in self.state.data_containers:        
      assert isinstance(a, TVArraySeries), type(a)
      r = TVDepsResource(a)
      aname = (a.name()[0], None)
      assert aname not in self.resources, namestr(aname)
      self.resources[aname] = r

    outputs = set()
    output_names = set()
          
    for tasknr, task in self.tasks.items():
      a = task.tvtask()
      for arg, argname in zip(a._inputs, a._inputs_names):
        if isinstance(arg, TVArray):
          assert arg.name() in self.data, namestr(arg.name())
          task.data.append(arg.name())
      if a._accumulate:
        for arg in a._outputs:
          prev = arg._previous()
          assert prev.name() in self.data, namestr(prev)
          task.data.append(prev.name())
      for arg, argname in zip(a._outputs, range(1, len(a._outputs)+1)):
        assert isinstance(arg, TVArray), (str(a), arg)
        assert arg.name() in self.data, namestr(arg.name())
        assert arg.name() not in output_names, namestr(arg.name())
        assert id(arg) not in outputs, namestr(arg.name())
        outputs.add(id(arg))
        output_names.add(arg.name())
        task.resources.append(arg.name()[:2])
        self.data[arg.name()].task = tasknr
    r = TVDepsResource(None)
    self.resources["@return"] = r
    for v in self.state.return_:
      mode, index = v
      assert mode == "data"
      a = self.state.data[index]
      if a._parent is not None and isinstance(a._parent, TVArraySeries):
        a = a._parent._children[-1]        
      assert a.name() in self.data, namestr(a)
      r.data.append(a.name())
            
    # 2: remove unassigned data
    for name0 in list(self.data.keys()):
      ok = False
      for name in list(self.data.keys()):
        if name[:2] != name0[:2]: continue
        if self.data[name].task is not None:
          ok = True
        if not ok:
          for r in self.resources.values():
            rdata = [d for d in r.data]
            if name in rdata:
              ok = True
              break
        if not ok:
          for r in self.tasks.values():
            rdata = [d for d in r.data]
            if name in rdata:
              ok = True
              break
        if ok: break    
      if not ok:
        self.data.pop(name0)
    
    # 3: check that every data has dtype and shape
    for d in self.data:
      r = self.data[d]
      tvdata = r.tvdata()
      if tvdata is None: continue
      assert tvdata.dtype is not None, d
      assert tvdata.shape is not None, d
    
    # 4: check consistency within containers
    for c in self.state.data_containers:
      assert isinstance(c, TVArraySeries)
      firstchild = c._children[0]
      for child in c._children[1:]:
        assert child.dtype == firstchild.dtype, (namestr(c.name()), child.dtype, firstchild.dtype)
        assert len(child.shape) == len(firstchild.shape), (namestr(c.name()), child.shape, firstchild.shape)
                                
    # 5: rewrite the depsgraph based on handlers
    for rname,r in list(self.resources.items()):      
      if rname not in self.resources: continue #if we were deleted by a previous handler
      if self.resources[rname] is not r: continue #if we were replaced by a previous handler
      tvdata = r.tvdata()
      if tvdata is None: continue
      if tvdata._parent is not None: continue
      assert isinstance(tvdata, TVArraySeries), r
      exceptions = []
      for handler in _handlers:
        try:
          handler(self, rname)
          break
        except:
          raise ###
          e = traceback.format_exc()
          exceptions.append(e)
      else:
        for e, handler in zip(exceptions, _handlers):
          print handler.__name__
          print e
          print
        raise Exception("All handlers failed for resource '%s'" % rname)
        

    # 6: add dependencies for shards, remove writeshards that point to eliminated data:
    for dname,d in list(self.data.items()):
      dname2 = dname[:2]      
      if isinstance(d.tvdata(), TVWriteShard):  
        if dname2 not in self.data:
          self.data.pop(dname)
          continue
        assert dname2 in self.data, dname2
        resname = (namestr(dname2)+"::join",None)
        if resname not in self.resources:
          self.resources[resname] = TVDepsResource(None)
          r = self.data[dname2].resource  
          if r is None:
            self.data[dname2].resource = resname
          else:
            self.resources[r].resources.append(resname)
          if d.tvdata()._parent._cache_active:
            cache = TVAllocator("cache", tvname=dname2)
            cachename = (namestr(dname2)+"::cache", None)
            self.allocators[cachename] = cache
            self.resources[resname].allocators.append(cachename)          
        self.resources[resname].data.append(dname)
      elif isinstance(d.tvdata(), TVReadShard):
        assert d.resource is None, (dname, d.resource)
        resname = (namestr(dname2)+"::split",None)
        if resname not in self.resources:
          self.resources[resname] = TVDepsResource(None)
          self.resources[resname].data.append(dname2)        
        d.resource = resname

    # 7: re-check that every data has dtype and shape (null data not allowed)
    for d in self.data:
      r = self.data[d]
      assert isinstance(r, TVDepsData), d
      tvdata = r.tvdata()
      assert tvdata is not None, d
      assert tvdata.dtype is not None, d
      assert tvdata.shape is not None, d

    # 8: calculate signatures for every task
    for task in self.tasks.values():
      t = task.tvtask()
      t.calc_signature()
          
    # 9: build hashtree for non cache-active data (such as those computed on the GPU)
    for d in self.data:
      r = self.data[d]      
      data = r.tvdata()
      if isinstance(data, TVReadShard): continue
      if not data._cache_active and (not hasattr(data, "_hash") or data._hash is None):
        if r.hashtree is None:
          r.hashtree = make_hashtree(d, self)
    
    
  def prune(self):
    """
    Eliminates all tasks with a cachekey, and everything that does not depend on return
    """
    #TODO: will run out of stack space for large depsgraphs
    tracing_data = set()
    tracing_resources = set()
    tracing_tasks = set()    
    dep_data = set()
    dep_resources = set()
    dep_tasks = set()
    dep_allocators = set()
    def trace_resource(name):
      if name in dep_resources: 
        return
      r = self.resources[name]
      tracing_resources.add(name)  
      for d in r.resources:
        if d in tracing_resources:
          raise ValueError("Resource '%s' and resource '%s' are in a cyclical dependency" % (name, d))
        trace_resource(d)
      for d in r.data:
        if d in tracing_resources:
          raise ValueError("Resource '%s' and data '%s' are in a cyclical dependency" % (name, d))        
        trace_data(d)
      for d in r.allocators:
        ok = True
        if self.cachedir is None:
          alloc = self.allocators[d]
          if alloc.alloc_command == "cache":
            ok = False
        if ok:    
          dep_allocators.add(d)
      dep_resources.add(name)   
      tracing_resources.remove(name)
    def trace_data(name):
      if name in dep_data: 
        return      
      r = self.data[name]
      if r.tvdata() is not None and r.tvdata()._cache.cachekey:
        dep_data.add(name) 
        for t in self.tasks:
          task = self.tasks[t]
          if name in task.data:
            task._cachekeys.append(r.tvdata())
        return
      tracing_data.add(name)
      d = r.task
      if d is not None:
        if d in tracing_tasks:
          raise ValueError("Data '%s' and task %d are in a cyclical dependency" % (namestr(name), d))        
        trace_task(d)              
      d = r.resource
      if d is not None:
        if d in tracing_resources:
          raise ValueError("Data '%s' and resource '%s' are in a cyclical dependency" % (namestr(name), d))        
        trace_resource(d)                      
      dep_data.add(name) 
      tracing_data.remove(name)
    def trace_task(name):
      if name in dep_tasks: 
        return
      r = self.tasks[name]
      tracing_tasks.add(name)  
      for d in r.resources:
        if d in tracing_resources:
          raise ValueError("Task %d and resource '%s' are in a cyclical dependency" % (namestr(name), d))
        trace_resource(d)
      for d in r.data:
        if d in tracing_resources:
          raise ValueError("Task %d and data '%s' are in a cyclical dependency" % (namestr(name), d))        
        trace_data(d)
      dep_tasks.add(name)
      tracing_tasks.remove(name)
      
    trace_resource("@return")
    for d in list(self.data.keys()):
      if d not in dep_data:
        self.data.pop(d)
    for d in list(self.resources.keys()):
      if d == "@return": continue
      if d not in dep_resources:
        self.resources.pop(d)
    for d in list(self.tasks.keys()):
      if d not in dep_tasks:
        self.tasks.pop(d)
    for d in list(self.allocators.keys()):
      if d not in dep_allocators:
        self.allocators.pop(d)

    retdep = self.resources["@return"].data
    progress = True
    while progress:
      progress = False
      for d in list(self.data.keys()):
        r = self.data[d]
        r.cleanup(self)      
        if d in retdep: continue
        if r.empty():
          progress = True
          self.data.pop(d)
      for d in list(self.resources.keys()):
        if d == "@return": continue
        r = self.resources[d]
        r.cleanup(self)
        if r.empty():
          progress = True 
          for a in r.allocators:
            self.allocators.pop(a)
          self.resources.pop(d)
      for d in list(self.tasks.keys()):
        r = self.tasks[d]
        r.cleanup(self)
        if r.empty():
          progress = True                    
          self.tasks.pop(d)      
    #print "PRUNE", len(self.tasks)
    
  def __str__(self):
    s = "@return:\n"
    namestrs = [namestr(dd) for dd in self.resources["@return"].data]
    s += "  " + " ".join(namestrs) + "\n"
    s += "Tasks:\n"
    for tnr,t in self.tasks.items():
      s += "  " + str(tnr) + ":\n"
      s += str(t)
    s += "Data:\n"
    for dname in sorted(list(self.data.keys())):
      d = self.data[dname]
      s += "  " + namestr(dname) + ":\n"
      s += str(d)      
    s += "Resources:\n"
    for rname in sorted(list(self.resources.keys())):
      r = self.resources[rname]
      if rname == "@return": continue
      s += "  " + namestr(rname) + ":\n"
      s += str(r)      
    s += "Allocators:\n"
    for aname in sorted(list(self.allocators.keys())):
      a = self.allocators[aname]
      s += "  " + namestr(aname) + ":\n"
      s += str(a)      
    return s
  
  def clone(self):
    newstate = self.state.clone()
    clone = TVDepsGraph(newstate, clone=True)    
    for k,v in self.tasks.items():
      d = v.dict(self.state)
      clone.tasks[k] = TVDepsTask.fromdict(d, clone.state)
    for k,v in self.data.items():
      d = v.dict(self.state)
      clone.data[k] = TVDepsData.fromdict(d, clone.state)
    for k,v in self.resources.items():
      d = v.dict(self.state)
      clone.resources[k] = TVDepsResource.fromdict(d, clone.state)
    for k,v in self.allocators.items():
      clone.allocators[k] = v.clone()
    for a in newstate.data:
      a._cache.graph = clone
    if self.cachedir is not None:
      clone.cachedir = self.cachedir
      clone.cache_indices = self.cache_indices
    return clone
    
  def define_cache(self, cachedir):
    from .TVCache import CacheIndex, cache_indices as global_cache_indices
    if cachedir is None: 
      self.cachedir = None
      self.cache_indices = None
      return
    self.cachedir = cachedir
    if not os.path.exists(self.cachedir):
      os.mkdir(self.cachedir)
    for d in self.data.values():
      data = d.tvdata()
      if data._cache is not None:
        data._cache.searched = False      
    cache_indices = {}
    for action in ("::shard","::join") + tuple(self.state.actions):
      k = (self.cachedir, action)
      if k not in global_cache_indices:
        dir = k[0] + os.sep + k[1]
        c = CacheIndex(dir)
        global_cache_indices[k] = c
      else:  
        c = global_cache_indices[k]
      cache_indices[action] = c
    self.cache_indices = cache_indices 
  
  def search_cache_hashes(self,allocator):
    from .TVArray import TVBaseShard, TVReadShard, TVWriteShard
    if self.cachedir is None: return
    global progress
          
    done = set()
    progress = True
    def search_cache_hash(dname):
      global progress
      if dname in done: return
      d = self.data[dname]
      data = d.tvdata()
      shards, shards2 = [], []
      if isinstance(data, TVBaseShard):
        dname2 = dname[:2]
        if isinstance(data, TVReadShard):          
          if dname2 not in done: return
          if data._parent._cache_active:
            data.search_cache_hash()
          done.add(dname)
          progress = True
          return      
      elif data.get_hash() is not None: 
        #print "HASH!", dname
        done.add(dname)
        progress = True
        return      
      else:        
        for sh in self.data:
          if len(sh) < 3: continue
          if sh[:2] != dname: continue
          shards.append(sh)
        shards2 = [sh for sh in shards if isinstance(self.data[sh].tvdata(), TVWriteShard)]  
        if d.task is None and not len(shards2):
          #print "UNUSED!", dname
          done.add(dname)
          progress = True
          return      
      if len(shards2): 
        #print "SEARCHIN?", dname
        deps = shards2
      else:  
        task = self.tasks[d.task]
        #print "SEARCHING", dname, str(task).split("\n")[0]
        deps = task.data
      if not d.tvdata()._cache_active:          
        ok = True
        for dep in deps:
          search_cache_hash(dep)
          if dep not in done:
            ok = False
            break
        if ok: 
          done.add(dname)  
          progress = True
        return
      for dep in deps:
        dd = self.data[dep]
        ddd = dd.tvdata()
        if isinstance(ddd, TVBaseShard):
          if dep not in done: return
        elif ddd.get_hash() is None:
          ok = False
          if not ddd._cache_active: 
            if dd.hashtree.bindable(self):
              ok = True
          if not ok:  
            return            
      #print "SEARCH", dname
      if len(shards2):
        data._cache.search_join()        
      else:
        data._cache.search()        
      if data._cache.arrayhash is not None:        
        #print "HASH2!", dname
        data.set_hash(data._cache.arrayhash)
        done.add(dname)
        progress = True
      #else:
      #  print "MISSED", task
    while progress:
      progress = False
      for d in self.data:
        search_cache_hash(d)
          
  def evaluate(self, cls_evaluator, args, kwargs):
    from . import get_mode, set_mode
    #TODO: check that all previous stages are complete
    assert not self._evaluated 
    self.state.bind_args(args, kwargs)
    
    evaluator = cls_evaluator(self)          
    self.data_unpruned = self.data.copy()
    
    oldmode = get_mode()
    set_mode("search_cache")
    self.search_cache_hashes(evaluator.allocator)    
    
    #prune depsgraph based on cachekeys
    self.prune()
    log.set_plan(str(self))
    
    try:      
      set_mode("evaluate")
      evaluator.init()
      evaluator.evaluate_resource("@return")
      evaluator.run()
      r = self.resources["@return"]
      ret = []
      for d in r.data:       
        dat = self.data[d]
        ret.append(dat.tvdata().get_data())
    finally:
      set_mode(oldmode)
      for d in self.data.values():
        dd = d.tvdata()
        if dd is not None:
          try:
            dd.set_data(None)      
          except AttributeError:
            pass
    self._evaluated = True  
    self.cache_indices = None
    ret = tuple(ret)
    if len(ret) == 1: 
      ret = ret[0]
    return ret  

_handlers = []
def add_handler(handler):
  """
  Adds a handler that can rewrite the depsgraph (especially resource dependencies)
  """
  from .TVArray import TVContainer
  assert callable(handler)
  _handlers.insert(0, handler)
  

def defaulthandler(graph, aname):
  """
  Default handler
  """
  r = graph.resources[aname]
  a = r.tvdata()

  ### 
  # Add a data dependency for a[n] on everything that depends on a[n-1] (except a[n] itself), including read shards!
  # Add a resource dependency for a[n] on a
  for cnr in range(a._maxchildindex):
    child = a._children[cnr]
    if child.name() == aname: continue
    rr = graph.resources[child.name()]
    rr.resources.append(aname)
    if cnr == 0: continue
    child0 = a._children[cnr-1]
    child0name0 = child0.name()
    child0names = [child0name0]
    for dname, d in graph.data.items():
      if len(dname) != 3: continue
      if dname[:2] != child0name0: continue
      if d.tvdata() is not None and isinstance(d.tvdata(), TVReadShard):
        child0names.append(dname)
    for taskname, task in graph.tasks.items():            
      if any([child0name in task.data for child0name in child0names]):
        for dname, d in graph.data.items():          
          if d.task == taskname and dname != child.name() and dname not in rr.data:
            rr.data.append(dname)
  if a._data is None:
    #Allocate data 
    gpu = ""
    if a.gpu: gpu = "gpu_"
    firstshape = a._children[0].shape
    for child in a._children[1:]:
      if child.shape != firstshape:
        same_shape = False
        break
    else:
      same_shape = True
    if same_shape:  
      #1. All data has the same shape, simple allocation 
      assert aname not in graph.allocators, aname
      alloc = TVAllocator(gpu+"alloc", tvname=aname, shape=a.shape, dtype=a.dtype)
      graph.allocators[aname] = alloc
      r.allocators.append(aname)
    else:
      #2. Different data shapes, malloc
      nbytes = 0
      for child in a._children:
        nbytes0 = functools.reduce(operator.__mul__, child.shape) * np.dtype(child.dtype).itemsize
        nbytes = max(nbytes, nbytes0)
      malloc = TVAllocator(gpu+"malloc", tvname=aname, nbytes=nbytes)
      graph.allocators[aname] = malloc
      r.allocators.append(aname)

          
from .TVArray import TVArray, TVArraySeries, TVBase, TVReadShard, TVWriteShard
from .TVAllocator import TVAllocator
add_handler(defaulthandler)