import numpy, sys, functools
from .TVCache import TVHash
from . import namestr, ViewError

class TVAllocator(object):  
  def __init__(self, alloc_command, **kwargs):
    self.alloc_command = alloc_command
    self._attrs = set()
    for k in kwargs:
      setattr(self, k, kwargs[k])      
  def __str__(self):
    s = "    " + self.alloc_command + "("
    first = True
    for k in sorted(list(self._attrs)):
      v = getattr(self, k)
      if v is None: continue
      if not first: s += ","
      first = False
      vv = str(v)
      if k == "dtype":
        try:
          vv = v.__name__
        except AttributeError:  
          vv = str(v) 
      elif k == "tvname" or k == "mallocname":
        vv = namestr(v)
      s += " " + k + " = " + vv
    s += ")\n"
    return s
  def dict(self):
    ret = {"alloc_command": self.alloc_command}
    for k in self._attrs:
      v = getattr(self, k)
      if v is not None:
        ret[k] = v
    return ret 
  def clone(self):
    clone = TVAllocator(self.alloc_command)
    for k in self._attrs:
      v = getattr(self, k)
      if v is not None:
        setattr(clone, k, v)
    return clone 
  def __setattr__(self, attr, value):
    if attr not in ("alloc_command", "_attrs"):
      self._attrs.add(attr)
    object.__setattr__(self, attr, value)
  def evaluate(self, allocator):
    assert isinstance(allocator, Allocator)
    kwargs = {}
    for k in list(self._attrs):
      kwargs[k] = getattr(self, k)
    allocator.command(self.alloc_command, kwargs)

def gpu_alloc(arr, name, dtype, shape):
  #print "GPU_ALLOC", name  
  import pycuda.gpuarray
  assert arr is not None, name
  m = pycuda.gpuarray.zeros(dtype=dtype, shape=shape)
  arr.set_data(m)

def gpu_malloc(arr, name, nbytes):
  #print "GPU_MALLOC", name  
  import pycuda.gpuarray
  assert arr is not None, name
  m = pycuda.gpuarray.zeros(nbytes, dtype=numpy.dtype('bool'))
  arr.set_data(m)

def alloc(arr, name, dtype, shape, memallocator = None):   
  #print "ALLOC", name  
  assert arr is not None, name
  if memallocator is not None:
    memallocate = memallocator.allocate
  else:
    memallocate = numpy.zeros
  data = memallocate(dtype=dtype, shape=shape)
  arr.set_data(data)

def malloc(arr, name, nbytes, memallocator = None): 
  #print "MALLOC", name  
  assert arr is not None, name
  if memallocator is not None:
    memallocate = memallocator.allocate
  else:
    memallocate = numpy.zeros
  data = memallocate(shape=(nbytes,),  dtype=numpy.dtype('bool'))
  arr.set_data(data)    
      
class Allocator(object):
  def __init__(self, graph, evaluator):
    self.graph = graph
    self.evaluator = evaluator
  def command(self, command, kwargs):
    if command == "alloc":
      return self.alloc(memallocator=self.evaluator.memallocator, **kwargs)
    elif command == "gpu_alloc":
      return self.gpu_alloc(**kwargs)    
    elif command == "malloc":
      return self.malloc(memallocator=self.evaluator.memallocator, **kwargs)        
    elif command == "gpu_malloc":
      return self.gpu_malloc(**kwargs)        
    elif command == "cache":
      return self.cache(**kwargs)        
    else:
      raise ValueError(command)
  def gpu_alloc(self, tvname, dtype, shape):
    data = self.graph.resources[tvname]
    arr = data.tvdata()
    gpu_alloc(arr, tvname, dtype, shape)        
  def alloc(self, tvname, dtype, shape, memallocator = None):
    data = self.graph.resources[tvname]
    arr = data.tvdata()
    alloc(arr, tvname, dtype, shape, memallocator)    
  def gpu_malloc(self, tvname, nbytes):
    data = self.graph.resources[tvname]
    arr = data.tvdata()
    gpu_malloc(arr, tvname, nbytes)    
  def malloc(self, tvname, nbytes, memallocator = None):
    data = self.graph.resources[tvname]
    arr = data.tvdata()
    malloc(arr, tvname, nbytes, memallocator)    
  def cache(self, tvname):
    o = self.graph.data[tvname].tvdata()
    if o._cache_active:
      #print "CACHE-SAVE2", o.name()
      arrayhash = o._cache.save_join(o.get_data(), o._cache_value)[1]
      outputhash = arrayhash
      if o.get_hash() is not None and o.get_hash() != outputhash: 
        print >> sys.stderr, "WARNING: Cache corruption detected for %s::join" % namestr(tvname)
      o.set_hash(outputhash)
    
def args_check_shapes(args):  
  for argname, arg in args.items():
    if not isinstance(arg, TVArray): continue
    assert arg.dtype is not None, argname
    assert arg.shape is not None, argname

def args_allocate(args):
  from . import get_mode
  assert get_mode() == "direct", get_mode()
  for argname, arg in args.items():
    if isinstance(arg, TVBaseShard):
      arg = arg._parent
    if not isinstance(arg, TVArray): continue
    assert arg.dtype is not None, argname
    assert arg.shape is not None, argname
    try:
      data = arg.get_data()
    except ViewError:
      #TODO: mechanism to ensure that arg is write-only (i.e. non-accumulating output)
      print "Warning: Need to delete and enlarge array '%s', hopefully it is write-only" % namestr(arg)
      data = None
    if data is None or arg.shape != data.shape: 
      if arg.gpu:
        gpu_alloc(arg, argname, arg.dtype, arg.shape)
      else:
        alloc(arg, argname, arg.dtype, arg.shape)

from .TVArray import TVArray, TVBaseShard
    