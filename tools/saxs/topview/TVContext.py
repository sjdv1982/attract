from . import namestr
import os, copy, inspect, numpy as np

class TVState(object):
  def __init__(self, func):
    self.func = func
    self.tasks = []
    self.data_containers = []
    self.data = []
    self.datadict = {}
    self.actions = set()
    self.return_ = None
    self.argbinding = None
    self._id_tasks = []
    self._id_data = []
    self._id_data_containers = []
  def find_task(self, task):
    return self._id_tasks.index(id(task))
  def find_data(self, data):
    return self._id_data.index(id(data))
  def find_data_container(self, data_container):
    return self._id_data_containers.index(id(data_container))
  def bind_template(self, register, args, kwargs, return_):
    """
    Populates tasks, data/containers and actions by providing the context register that was filled by an "unroll" call to self.func,
     using 'return_ = self.func(*args, **kwargs)'
    Binds all TVArrays built from input args/kwargs to the appropriate arg
    """
    from .TVAction import TVTask
    from .TVArray import TVArray, TVBaseShard
    self.tasks = []
    self.data_containers = []
    self.data = []
    self.datadict = {}
    self.actions = {}
    self._id_tasks = []
    self._id_data = []
    self._id_data_containers = []
    for a in register:
      if isinstance(a, TVTask):
        self.tasks.append(a)
        self._id_tasks.append(id(a))
      elif isinstance(a, TVArray):
        self.data.append(a)
        self._id_data.append(id(a))
        assert a.name() not in self.datadict, (namestr(a.name()), a, self.datadict[a.name()], a._pindex, self.datadict[a.name()]._pindex)
        self.datadict[a.name()] = a
    def add_container(a):
      if isinstance(a, TVBaseShard): return
      p = a._parent
      if p is None: return
      if id(p) in self._id_data_containers: return
      self.data_containers.append(p)
      self._id_data_containers.append(id(p))
      add_container(p)
    for a in self.data:
      add_container(a)
    for a in self.tasks:
      self.actions[a.action.name()] = a.action
    
    self.argbinding = {} #data index => (function arg name, arg access index, arg length)
    callargs = inspect.getcallargs(self.func, *args, **kwargs)
    for anr, a in enumerate(self.data):      
      if isinstance(a, TVBaseShard): continue
      data = a.get_data() 
      if data is None: continue      
      ok = False
      if a._accessor is not None:        
        accessor, index = a._accessor
        acdata = accessor._data
        hash = None
        for k,v in callargs.items():
          if isinstance(v, TVHash) and v is a._hash:   
            hash = k          
          elif v is acdata:
            self.argbinding[anr] = (k, index, len(acdata), hash)
            ok = True
            break        
      else:
        hash = None        
        arg = None
        for k,v in callargs.items():
          if isinstance(v, TVHash) and v is a._hash:   
            hash = k
          elif v is data:
            arg = k
            ok = True
        if ok:
          self.argbinding[anr] = (arg,  None, None, hash)
      if not ok:
        raise ValueError("Cannot bind the array data of TVArray '%s' to any function argument" % namestr(a.name()))
  
    #bind return
    ret = []  
    for argnr, arg in enumerate(return_):
      if isinstance(arg, TVArray):        
        try:
          anr = self.find_data(arg)
        except ValueError:
          raise ValueError("Cannot bind return argument %d to any TVArray from the function" % (argnr+1))        
        ret.append(("data", anr))          
      #TODO: allow containers as well (check in self.data_containers)
    self.return_ = tuple(ret)    
  
    #determine list of arguments that don't bind (and therefore may determine control flow)
    bindargs = [v[0] for v in self.argbinding.values()] + [v[3] for v in self.argbinding.values() if v[3] is not None]  
    nonbind = sorted([k for k in callargs if k not in bindargs])
    return nonbind
  
  def bind_actions(self, actiondict):
    for a in self.tasks:
      action = actiondict[task._actionname]
      assert action.name() == task._actionname #an entry 'foo' in the action dict must have 'foo' as its name()
      task.action = action
      task._actionname = None
  def bind_args(self, args, kwargs):
    callargs = inspect.getcallargs(self.func, *args, **kwargs)
    for anr, k in self.argbinding.items():
      data = self.data[anr]
      assert data is not None, data.tvname
      bindarg = callargs[k[0]]
      index, length, hash = k[1], k[2], k[3]
      if index is not None and length is not None:
        assert len(bindarg) == length, (k[0], len(bindarg), length)
        bindarg = bindarg[index]      
      assert data == bindarg
      if data._gpu:
        import pycuda.gpuarray
        assert isinstance(bindarg, pycuda.gpuarray.GPUArray), k
        assert hash is not None, k
        data.set_hash(hash)
      else:  
        assert isinstance(bindarg, np.ndarray), k      
        data.set_hash(TVHash(bindarg))
      data.set_data(bindarg)
  def dict(self):
    assert self.return_ is not None
    assert self.argbinding is not None    
    tasks = [t.dict(self) for t in self.tasks]
    data = [d.dict(self) for d in self.data]
    data_containers = [d.dict(self) for d in self.data_containers]
    data_container_types = [d.__class__.__name__ for d in self.data_containers]
    actions = list(self.actions.keys())
    return {
      "tasks": tasks,
      "data": data,
      "data_containers": data_containers,
      "data_container_types": data_container_types,
      "actions": actions,
      "return": self.return_,
      "argbinding": self.argbinding
    }
  def set_actions(self, actions):
    from .TVAction import TVAction
    assert isinstance(actions, dict)
    for k,v in actions.items():
      assert isinstance(k, str), k
      assert isinstance(v, TVAction), v
    self.actions = actions
  def set_dict(self, dic):
    from .TVArray import TVArray, TVArraySeries, TVReadShard, TVWriteShard
    from .TVAction import TVTask
    assert self.func is not None
    assert self.actions is not None
    for action in dic["actions"]:
      assert action in self.actions, action
    types = {
      "TVArraySeries": TVArraySeries,
    }
    self.tasks = []
    self.data = []
    self.data_containers = []
    self.datadict = {}
    self._id_data = []
    self._id_data_containers = []
    self._id_tasks = []
    for typename, ddic in zip(dic["data_container_types"], dic["data_containers"]):
      assert isinstance(ddic, dict), ddic
      typ = types[typename]
      container = typ.fromdict(ddic, self)
      self.data_containers.append(container)
      self._id_data_containers.append(id(container))
    
    self.data = [None for ddic in dic["data"]]
    self._id_data = [None for ddic in dic["data"]]
    for ddnr, ddic in enumerate(dic["data"]):
      assert isinstance(ddic, dict), ddic
      if "offset" in ddic: continue
      data = TVArray.fromdict(ddic, self)
      self.data[ddnr] = data
      self.datadict[data.name()] = data
      self._id_data[ddnr] = id(data)
    for ddnr, ddic in enumerate(dic["data"]):
      assert isinstance(ddic, dict), ddic
      if "offset" not in ddic: continue
      if ddic["write"]:
        data = TVWriteShard.fromdict(ddic, self)
      else:
        data = TVReadShard.fromdict(ddic, self)
      self.data[ddnr] = data
      self.datadict[data.name()] = data
      self._id_data[ddnr] = id(data)
    for ddic in dic["tasks"]:
      assert isinstance(ddic, dict), ddic
      task = TVTask.fromdict(ddic, self)
      self.tasks.append(task)
      self._id_tasks.append(id(task))
    self.argbinding = dic["argbinding"]
    self.return_ = dic["return"]  
  def clone(self):
    from .TVArray import TVArray, TVContainer
    clone = TVState(self.func)
    clone.set_actions(self.actions)
    
    d = self.dict()
    clone.set_dict(d)    
    return clone
    
    
class TVContext(object):
  def __init__(self, func, cachedir=None, evaluator=None):
    self.func = func
    self._graphs = {}
    self._nonbind = None
    if cachedir is not None:
      cachedir = os.path.abspath(cachedir)
    self._cachedir = cachedir
    
    from .Evaluator import Evaluator, SimpleEvaluator
    if evaluator is not None:
      assert issubclass(evaluator, Evaluator)
    else:
      evaluator = SimpleEvaluator
    self.cls_evaluator = evaluator
    
  def calc_signature(self, *args, **kwargs):  
    from .TVCache import TVHash
    assert self._nonbind is not None
    callargs = inspect.getcallargs(self.func, *args, **kwargs)
    hashes = []
    for a in self._nonbind:
      arg = callargs[a]
      if isinstance(a, np.ndarray):
        arg = TVHash(arg)
      hashes.append(hash(arg))
    return tuple(hashes)  
    
  def unroll(self, *args, **kwargs):
    from . import get_mode, set_mode
    from .TVArray import TVData    
    oldmode = get_mode()
    set_mode("unroll")
    _push(self)
    self._reg  = []
    ret = self.func(*args, **kwargs)
    assert isinstance(ret, tuple) or isinstance(ret, TVData), ret
    if isinstance(ret, TVData): ret = (ret,)
    state = TVState(self.func)
    nonbind = state.bind_template(self._reg, args, kwargs, ret)
    c = _pop()
    assert c is self
    set_mode(oldmode)
    graph = TVDepsGraph(state)
    if self._nonbind is None:
      self._nonbind = nonbind
    else:
      assert nonbind == self._nonbind, (nonbind, self._nonbind) 
    return graph
  
  def get_graph(self, *args, **kwargs):
    if self._nonbind is None:
      graph = self.unroll(*args, **kwargs)
      signature = self.calc_signature(*args, **kwargs)      
      self._graphs[signature] = graph
    else:
      signature = self.calc_signature(*args, **kwargs)
      graph = self._graphs.get(signature, None)
      if graph is None:
        graph = self.unroll(*args, **kwargs)
        self._graphs[signature] = graph
    if graph.cachedir != self._cachedir:
      graph.define_cache(self._cachedir)    
    return graph
  
  def __call__(self, *args, **kwargs):    
    graph0 = self.get_graph(*args, **kwargs)
    graph = graph0.clone()
    result = graph.evaluate(self.cls_evaluator, args, kwargs)
    graph = None
    import gc
    gc.collect()
    return result
    
    
_stack = []
def _push(context):
  _stack.append(context)

def _pop():
  assert len(_stack)
  return _stack.pop()
  
def register(obj):
  from .TVAction import TVTask
  from .TVArray import TVContainer, TVArray
  assert isinstance(obj, TVTask) or isinstance(obj, TVContainer) or isinstance(obj, TVArray), obj
  assert len(_stack)
  _stack[-1]._reg.append(obj)

from .TVCache import TVHash
from .TVDepsGraph import TVDepsGraph