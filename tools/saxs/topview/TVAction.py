from . import namestr, TVBase, TVData, get_mode 
from .TVContext import register
from .TVArray import TVArrayBase, TVArray, TVBaseShard, TVReadShard, TVWriteShard
from .TVAllocator import args_check_shapes, args_allocate
from collections import namedtuple, OrderedDict
import numpy as np
import sys, functools, copy
has_spyder =False
#try:
  #import spyder, Spyder
  #has_spyder = True
#except ImportError:
  #has_spyder = False  


def valid_argument(a):
  if isinstance(a, TVArrayBase): return True
  if isinstance(a, TVBaseShard): return True
  if isinstance(a, (str, int, bool, float, np.generic)): return True
  if isinstance(a, tuple): 
    if not all([valid_argument(aa) for aa in a]):
      return False
    return True
  if has_spyder and isinstance(a, Spyder.Object):
    return True
  return False

def serialize_argument(a, namedtupledict):
  if isinstance(a, TVArray):
    return "TVArray", a.name()
  if isinstance(a, (str, int, bool, float)):
    return "builtin", a
  if isinstance(a, np.generic):
    return "binary", (type(a).__name__, a.tostring())
  if isinstance(a, tuple):
    d = tuple([serialize_argument(aa) for aa in a])
    if isinstance(a, namedtuple):
      cls = a.__class__
      clsname = cls.__name__
      if clsname in namedtupledict:
        cls2 = namedtupledict[clsname]
        assert cls._fields == cls2._fields, (clsname, cls._fields, cls2._fields) #namedtuples with the same name must have the same fields
      namedtupledict[clsname] = cls
      return "namedtuple", (clsname, a._fields, d)
    else:
      return "tuple", d
  if has_spyder and isinstance(a, Spyder.Object):
    try:
      d = a.dict()
    except AttributeError:
      d = str(a)
    return "Spyder", (a.typename(), d)
    
def deserialize_argument(a, state, namedtupledict):
  mode, data = a
  if mode == "TVArray":
    return state.datadict[data]
  if mode == "builtin":
    return data
  if mode == "binary":
    typ = numpy.getattr(data[0])
    return np.fromstring(data[1], dtype=typ)[0]
  if mode == "tuple":
    return tuple([deserialize_argument(dd, state, namedtupledict) for dd in data])
  if mode == "namedtuple":
    name, fields = data[0], data[1]
    if name not in namedtupledict:
      namedtupledict[name] = namedtuple(name, fields)
    cls = namedtupledict[name]  
    d = [deserialize_argument(dd, state, namedtupledict) for dd in data[2]]
    return cls(d)
  if mode == "Spyder":
    assert has_spyder
    cls = Spyder.load(data[0])
    try:
      d = cls.fromdict(data[1])
    except AttributeError:
      d = cls(data[1])
    return d  
  raise TypeError(mode)
  
def convert_argument(a):
  if isinstance(a, TVArray):
    data = a.get_data()
    assert data is not None 
    assert a.shape == data.shape, (a.shape, data.shape)
    assert a.dtype == data.dtype, (a.dtype, data.dtype)
    a = data
  elif isinstance(a, TVBase):
    raise TypeError(a.__class__)
  return a
 
class TVTask(TVBase):
  def __init__(self, action, input_args, input_kwargs):
    assert isinstance(action, TVAction)
    self.action = action
    self._input_args, self._input_kwargs = self._current_arguments(input_args, input_kwargs)        
    self._inputs = self._input_args + tuple([self._input_kwargs[k] for k in sorted(self._input_kwargs.keys())])
    for i in self._inputs:
      if isinstance(i, TVWriteShard):
        raise TypeError("Cannot use write shard '%s' as input", namestr(i.name()))
      elif isinstance(i, TVReadShard):
        assert i._parent._current() is i._parent, ( namestr(i._parent._current().name()), namestr(i._parent.name()), namestr(i.name()))
      elif isinstance(i, TVArray):
        assert i.is_assigned() or i.get_data() is not None, namestr(i.name())
    self._inputs_names = tuple(range(1,len(self._input_args)+1)) + tuple(sorted(self._input_kwargs.keys()))
    
  @staticmethod
  def _current_arguments(args, kwargs):    
    args2, kwargs2 = [], {}
    for a in args:
      assert valid_argument(a), (type(a), isinstance(a, TVArrayBase), isinstance(a, TVArray), a)
      if isinstance(a, TVArrayBase):
        a = a._current()
      elif isinstance(a, TVBaseShard):  
        pass
      else:
        a = copy.deepcopy(a)
      args2.append(a)  
    for name, a in kwargs.items():
      assert valid_argument(a), name
      if isinstance(a, TVArrayBase):
        a = a._current()
      elif isinstance(a, TVBaseShard):  
        pass
      else:
        a = copy.deepcopy(a)
      kwargs2[name] = a
    return tuple(args2), kwargs2 
  def __gt__(self, output):
    return self._assign(output, accumulate=False)
  def __rshift__(self, output):
    return self._assign(output, accumulate=True)  
  def _assign(self, output, accumulate):  
    if isinstance(output, TVData):
      output = (output,)
    else:
      assert isinstance(output, tuple), output
    
    for o in output:
      assert isinstance(o, TVData), o
      for i in self._input_args + tuple(self._input_kwargs.values()):
        if i is o: 
          raise ValueError("Action '%s': TVArray '%s' appears in both input and output" % (self.action.name(), namestr(i.name())))            
    output0 = []
    for o in output:
      if get_mode() != "direct":
        if accumulate and not o._assigned:
          raise ValueError("Cannot accumulate to unassigned array '%s'" % namestr(o.name()))
      oo = o._assign(accumulate=accumulate)
      assert isinstance(oo, TVData), oo
      output0.append(oo)
    output = tuple(output0)

    self._outputs = self._current_arguments(output, {})[0]
    self._accumulate = accumulate
    self.args = self.action.format_arguments(self._input_args, self._input_kwargs, self._outputs, self._accumulate)    
    assert isinstance(self.args, OrderedDict)
    self.action.shape(self.args)
    args_check_shapes(self.args)
     
    mode = get_mode()
    if mode == "direct":
      args_allocate(self.args)      
      #print >> sys.stderr, "EVALUATE", self
      self.action._evaluate_direct(self.args, self._outputs)
    elif mode  == "unroll":            
      register(self)
    else:
      raise ValueError(mode)
    
  def calc_signature(self):
    """
    The purpose of signatures is to store the hashes of parameters, allowing for optimization
    """
    sig = []    
    for arg in self.args.values():
      if isinstance(arg, TVArray):
        params = arg.params(hashing=True)
        try:
          h = hash(tuple(params.iteritems()))
        except TypeError:
          raise TypeError("Unhashable params in TVArray '%s': '%s'" % (namestr(arg), params))
      elif isinstance(arg, (str, int, bool, float, np.generic)): 
        h = arg
      else:
        h = hash(arg)        
      sig.append(h)
    self.signature = sig  
    
  def evaluate(self, cache):  
    #print "EVALUATE", self
    self.action._evaluate(self.args)
    if cache and get_mode() == "evaluate":
      for o in self._outputs:
        if o._cache_active:
          #print "CACHE-SAVE", o.name()
          arrayhash = o._cache.save(o.get_data(), o._cache_value)[1]
          outputhash = arrayhash
          if o.get_hash() is not None and o.get_hash() != outputhash: 
            print >> sys.stderr, "WARNING: Cache corruption detected for action '%s' producing '%s'" % (self.action.name(),namestr(o.name()))
          o.set_hash(outputhash)

  def __str__(self):
    n = self.action.name()
    s = n + "("
    first = True
    for anr, a in enumerate(self._input_args):
      if not first:
        s += ", "
      first = False  
      if isinstance(a, TVData):
        aname = namestr(a.name())
        if aname is None: raise TypeError(str(a))
        s += aname        
      else:
        s += str(a)
    for aname, a in self._input_kwargs.items():
      if not first:
        s += ", "
      first = False  
      s += "%s = " % aname
      if isinstance(a, TVData):
        s += namestr(a.name())
      else:
        s += str(a)
    token = ">"
    if self._accumulate: token = ">>"
    s += ") %s " % token      
    if len(self._outputs) != 1:
      s += "("
    for anr, a in enumerate(self._outputs):
      if anr > 0:
        s += ", "
      if isinstance(a, TVData):
        s += namestr(a.name())
      else:
        s += str(a)
    if len(self._outputs) != 1:
      s += ")"      
    return s
  def _dict(self, state, serialize, namedtupledict=None):
    d = {}
    def dict2(a):
      if serialize:
        return serialize_argument(tuple(aa), namedtupledict)
      else:
        ret = []
        for aa in a:
          if isinstance(aa, TVArray):
            assert isinstance(aa, TVBaseShard) or hasattr(aa, "_tvname"), aa
            v = True, state.find_data(aa)
          else:
            v = False, aa
          ret.append(v)
        return ret
    d["title"] = str(self)
    d["action"] = self.action.name()
    d["_input_args"] = dict2(self._input_args)
    d["_input_kwargs"] = dict(zip(self._input_kwargs.keys(), dict2(self._input_kwargs.values())))
    d["_outputs"] = dict2(self._outputs)
    d["_accumulate"] = self._accumulate
    assert self.signature is not None
    d["signature"] = self.signature
    return d
  def dict(self, state):
    return self._dict(state, serialize=False)
  def serialize(self, namedtupledict):
    return self._dict(state, serialize=True, namedtupledict=namedtupledict)
  
  @classmethod
  def _fromdict(cls, dic, state, deserialize, namedtupledict=None):
    self = object.__new__(cls)
    self.action = state.actions[dic["action"]]
    def fromdict2(a):
      if deserialize:
        return deserialize_argument(a, state, namedtupledict)
      else:
        ret = []
        for is_array, aa in a:
          if is_array:
            v = state.data[aa]
          else:
            v = aa
          ret.append(v)
        return tuple(ret)  
    self._input_args = fromdict2(dic["_input_args"])
    ik = dic["_input_kwargs"]
    input_kwarg_values = fromdict2(tuple(ik.values()))
    self._input_kwargs = dict(zip(ik.keys(), input_kwarg_values))
    self._inputs = self._input_args + tuple([self._input_kwargs[k] for k in sorted(self._input_kwargs.keys())])
    self._inputs_names = tuple(range(1,len(self._input_args)+1)) + tuple(sorted(self._input_kwargs.keys()))
    
    self._outputs = fromdict2(dic["_outputs"])
    self._accumulate = dic["_accumulate"]
    self.args = self.action.format_arguments(self._input_args, self._input_kwargs, self._outputs, self._accumulate)    
    args_check_shapes(self.args)
    self.signature = dic["signature"]    
    return self      
  
  @classmethod
  def fromdict(cls, dic, state):
    return cls._fromdict(dic, state, deserialize=False)
  @classmethod
  def deserialize(cls, dic, state, namedtupledict=None):
    if namedtupledict is None: namedtupledict = {}
    return cls._fromdict(dic, state, deserialize=True, namedtupledict=namedtupledict)
                    
class TVAction(TVBase):
  def run(self, args):
    raise NotImplementedError 
  def shape(self, args):
    raise NotImplementedError
  def validate(self, args):
    raise NotImplementedError  
  def format_arguments(self, input_args, input_kwargs, output, accumulate):
    raise NotImplementedError
   
  def _evaluate(self, args):
    args = self._convert(args) 
    self.run(args)
  
  def _evaluate_direct(self, args, outputs):
    return self._evaluate(args)
    
  @classmethod
  def _convert(cls, args):    
    args2 = OrderedDict()
    for aname, a in args.items():
      try:
        a = convert_argument(a)
      except Exception as e:
        raise TypeError, "Action '%s': Cannot convert argument '%s'" % (cls.__name__, aname) + "\n  " + sys.exc_info()[0].__name__ + " : " + str(sys.exc_info()[1]), sys.exc_info()[2]
      args2[aname] = a
    return args2
    
  def __new__(cls, *args, **kwargs):    
    self = object.__new__(cls)
    return TVTask(self, args, kwargs)
  
  @classmethod
  def name(cls):
    return cls.__name__
  
class GPUGrid(object):
  def __init__(self, grid, block):
    assert len(grid) <= 2
    for n in range(len(grid)): 
      assert isinstance(grid[n], int), n
    assert len(block) <= 3
    for n in range(len(block)): 
      assert isinstance(block[n], int), n      
    self.block = block
    self.grid = grid
    
class GPUPrepare(object):
  def __init__(self, arg_types, texrefs=[], shared_size=0):
    self.arg_types = arg_types
    self.texrefs = texrefs
    self.shared_size = shared_size
    
class GPUAction(TVAction): 
  _gpu_mod_cache = None
  _gpu_func_cache = None  
  source = None
  main = None
  stream = None
  def define_grid(self, args):
    raise NotImplementedError
  def prepare(self, funcname, args):
    raise NotImplementedError
  def get_source(self, args):
    assert self.source is not None
    return self.source
  def convert(self, args):
    import pycuda
    ret = OrderedDict()
    for aname, a in args.items():
      if isinstance(a, TVArray):
        assert a.gpu == True, aname
        assert a.get_data() is not None, aname
        a = a.get_data()
      elif isinstance(a, TVBase):
        raise TypeError((aname, a.__class__))
      if isinstance(a, pycuda.gpuarray.GPUArray):
        a = a.gpudata
      ret[aname] = a
    return ret  
  def kernel_args(self, args):
    return tuple(args.values())
      
  def __new__(cls, *args, **kwargs):    
    import pycuda.autoinit  
    if cls._gpu_func_cache is None: 
      cls._gpu_func_cache = {}
    if cls._gpu_mod_cache is None: 
      cls._gpu_mod_cache = {}    
    assert cls.main is not None and isinstance(cls.main, str)
    return super(GPUAction, cls).__new__(cls, *args, **kwargs)    
  
  def _evaluate(self, args):
    import pycuda
    from pycuda.autoinit import context
    from pycuda.compiler import SourceModule
    source = self.get_source(args)
    if not isinstance(source, (str, bytes)):
      raise TypeError("GPUAction '%s': self.get_source must return a str or bytes instance" % self.__class__.__name__)      
    mod = self._gpu_mod_cache.get(source, None)
    if mod is None:
      mod = SourceModule(source)
      self._gpu_mod_cache[source] = mod
    prepared_main = self._gpu_func_cache.get((source, self.main), None)  
    if prepared_main is None:
      main = mod.get_function(self.main)  
      gpu_prepare = self.prepare(main, args)
      if not isinstance(gpu_prepare, GPUPrepare):
        raise TypeError("GPUAction '%s': self.prepare must return a GPUPrepare instance" % self.__class__.__name__)      
      prep = main.prepare(gpu_prepare.arg_types, texrefs=gpu_prepare.texrefs)      
      prepared_main = prep, gpu_prepare.shared_size
      self._gpu_mod_cache[source, self.main] = prepared_main
    gpu_grid = self.define_grid(args)
    if not isinstance(gpu_grid, GPUGrid):
      raise TypeError("GPUAction '%s': self.define_grid must return a GPUGrid instance" % self.__class__.__name__)          
    conv_args = self.convert(args)
    if not isinstance(args, OrderedDict):
      raise TypeError("GPUAction '%s': self.convert must return a OrderedDict instance" % self.__class__.__name__)          
    kernel_args = self.kernel_args(conv_args)
    assert isinstance(kernel_args, tuple), kernel_args
    
    prep, shared_size = prepared_main
    if self.stream is None:      
      prep.prepared_call(gpu_grid.grid, gpu_grid.block,*kernel_args, shared_size=shared_size)
    else:
      prep.prepared_async_call(gpu_grid.grid, gpu_grid.block,self.stream,*kernel_args, shared_size=shared_size)
    
    self.post_evaluate(args)
    
  def post_evaluate(self, args):
    pass
  
  def set_stream(self, stream):
    from pycuda.driver import Stream
    assert stream is None or isinstance(stream, Stream)
    self.stream = stream
    
    
    
from .TVCache import TVHash 