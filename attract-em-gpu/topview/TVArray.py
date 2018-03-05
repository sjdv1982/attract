from . import namestr, TVBase, TVData, get_mode, ViewError
import numpy as np
from collections import OrderedDict
import copy

try:
  import pycuda.gpuarray
  has_pycuda = True
except ImportError:
  has_pycuda = False

class TVContainer(TVData):
  _parent = None  
  _index = None
  _tvname = None
  _gpu = False
  _gpucontext = None
  _attrs = ("dtype", "flags", "shape", "strides", "_gpu", "_gpucontext")
  def _set_parent(self, parent, index):
    assert parent._tvname == self._tvname
    self._parent = parent
    self._index = index
  def name(self):
    if self._parent is None or self._index is None: 
      return (self._tvname, None)
    else:
      return (self._tvname, str(self._index))
  def __eq__(self, other):
    types = (TVContainer, np.ndarray)
    if has_pycuda:
      types = types + (pycuda.gpuarray.GPUArray,)
    if not isinstance(other, types): return False
    for a in self._attrs:
      if a.startswith("_"): continue
      v1 = getattr(self, a)
      v2 = getattr(other, a)
      if a == "flags" and v1 is None: continue
      if a == "strides" and v1 is None: 
        v1 = self.get_strides()
      if a == "dtype":
        v1 = np.dtype(v1)
        v2 = np.dtype(v2)
      if v1 != v2:
        return False
    return True  
  def __ne__(self, other):
    return not self.__eq__(other)    
 
class TVArrayList(TVData):
  """
  Structure to bind a list of numpy arrays to TVArrays
  Only exists during in unroll and direct mode 
  Not maintained in depsgraph: stored as accessor instead
  """
  def __init__(self, tvname, data, hashes = None, **kwargs):
    assert "hash" not in kwargs, list(kwargs.keys())
    assert isinstance(tvname, str), tvname        
    self._data = data
    self._items = []
    for dnr in range(len(data)):
      subdata = data[dnr]
      hash = None
      if hashes is not None:
        hash = hashes[dnr]
      item = TVArray(tvname+"[%d]" % dnr, subdata, __accessor__= (self, dnr), hash=hash, **kwargs)
      self._items.append(item)
  def __getitem__(self, itemnr):
    return self._items[itemnr]
  def __len__(self):
    return len(self._items)
  def __iter__(self):
    return self._items.__iter__()

class TVArrayContainer(TVContainer):
  _maxchildindex = 0
  def __init__(self, tvname):
    assert isinstance(tvname, str), tvname
    self._tvname = tvname
    self._children = []
    if get_mode() == "unroll":
      register(self)    
  def __len__(self):
    assert self._maxchildindex == len(self._children)
    return self._maxchildindex
    
class TVArraySeries(TVArrayContainer):
  """
  Never explicitly created by the user, but always automatically created as a parent for a TVArray"
  """
  _cache_value = False
  _cache_active = True
  def __init__(self, array):
    TVArrayContainer.__init__(self, array._tvname)
    assert isinstance(array, TVArray)
    self._maxchildindex = 1
    self._children.append(array)  
  def clear(self):
    self._children[0].clear()
  def __getitem__(self, index):
    assert isinstance(index, int)
    assert index >= 0
    assert get_mode() not in ("direct", "unroll")
    return self._children[index]
  def _assign(self, accumulate):
    assert len(self._children)
    current = self._children[-1]
    child = TVArray(self._tvname, __parent__ = True)
    for attr in self._attrs:
      value = getattr(current, attr)      
      setattr(child, attr, value)          
    child._set_parent(self, self._maxchildindex)
    self._children.append(child)
    self._maxchildindex += 1    
    child._cache_value = self._cache_value
    child._cache_active = self._cache_active
    return child
  def __getattr__(self, attr):
    return getattr(self._children[0], attr)
  def dict(self, state):
     ret = {
       "_tvname": self._tvname,
       "_cache_value": self._cache_value,
       "_cache_active": self._cache_active,
       "_maxchildindex": self._maxchildindex,
     }  
     return ret
  @classmethod
  def fromdict(cls, dic, state):
    self = object.__new__(cls)
    TVArrayContainer.__init__(self, dic["_tvname"])
    self._cache_value = dic["_cache_value"]
    self._cache_active = dic["_cache_active"]
    self._maxchildindex = dic["_maxchildindex"]
    self._hash = None
    return self

class _MetaArray(type):
  def __instancecheck__(cls, other):
    if isinstance(other, TVBaseShard):
      return True
    return type.__instancecheck__(cls, other)
  
class TVArrayBase(TVContainer):
  __metaclass__ = _MetaArray

class TVArray(TVArrayBase):      
  _data = None
  _assigned = False
  _cache = None
  _hash = None
  _cache_value = False
  _cache_active = True
  _join_assigned = False
  _sharded = False
  def __init__(self, tvname, *args, **kwargs):
    from .TVCache import TVHash
    assert isinstance(tvname, str), tvname
    self._tvname = tvname
    assert len(args) in (0,1), args
    if len(args) == 0 and "data" in args:
      args = [kwargs.pop("data")]
    if len(args) == 1:
      a = args[0]
      gpu = False
      if "gpu" in kwargs:
        gpu = kwargs["gpu"]
      if gpu:
        assert has_pycuda
        if isinstance(a, np.ndarray):
          d = pycuda.gpuarray.to_gpu(a)
          if get_mode() != "direct":
            self._hash = TVHash(d)
        elif isinstance(a, pycuda.gpuarray.GPUArray):
          d = a
          if get_mode() != "direct":
            if "hash" not in kwargs:
              raise TypeError("TVArrays can be constructed from data on the GPU, but then you have to supply a TVHash hash argument, computed before it was sent to the GPU")
            hash = kwargs["hash"]
            assert isinstance(hash, TVHash)
            self._hash = hash
        else:
          raise TypeError
        self._data = d
      else:  
        assert isinstance(a, np.ndarray)
        self._data = a
      for p in self._attrs:
        pp = p
        if p[0] == "_": pp = p[1:]
        if pp != p or not hasattr(a, p):
          setattr(self, p, None)
        elif pp in kwargs:
          setattr(self, p, kwargs[pp])
        else:
          setattr(self, p, getattr(a, p))
      self._gpu = gpu    
    else:            
      for p in self._attrs:
        pp = p
        if p[0] == "_": pp = p[1:]        
        setattr(self, p, None)
        if pp in kwargs:
          setattr(self, p, kwargs[pp])
    if get_mode() == "unroll":
      register(self)
    if self._gpu:
      self._cache_active = False  
    if "__parent__" not in kwargs:
      self._parent = TVArraySeries(self)
      self._index = 0
      self._parent._cache_active = self._cache_active
    self._cache = TVCache(self) 
    self._accessor = kwargs.get("__accessor__", None)
  def get_strides(self):    
    if hasattr(self, "strides") and self.strides is not None:
      return self.strides
    data = self.get_data() 
    if data is not None:
      return data.strides
    assert self.flags is None or self.flags.c_contiguous
    itemsize = np.dtype(self.dtype).itemsize
    factor = 1
    strides = []
    for n in reversed(range(len(self.shape))):
      strides.append(factor*itemsize)
      factor *= self.shape[n]
    strides = tuple(reversed(strides))
    return strides
  def name(self):
    if len(self._parent._children) == 1: 
      return (self._tvname, None)
    else:
      return (self._tvname, str(self._index))
    
  def params(self, hashing=False):
    ret = OrderedDict()
    for p in self._attrs:
      data = self.get_data()
      if data is not None and hasattr(data, p):
        pp = getattr(data, p)      
      else:  
        pp = getattr(self, p)      
      ppp = None
      if pp is not None:        
        ppp = pp
        if p == "flags":          
          if not pp.c_contiguous:
            pass
          elif not hasattr(pp, "owndata") or not pp.owndata:
            pass
          elif not pp.writeable:
            pass
          elif not pp.aligned:
            pass
          elif pp.updateifcopy:  
            pass
          else:
            ppp = None
        if hashing and ppp is not None:
          ppp = str(ppp)
      if ppp is not None:   
        ret[p] = ppp
    return ret
    
  def __str__(self):
    ret = namestr(self.name()) + " = TVArray("    
    first=True
    for p in self._attrs:
      data = self.get_data()
      if data is not None and hasattr(data, p):
        pp = getattr(data, p)      
      else:  
        pp = getattr(self, p)      
      ppp = None
      if pp is not None:        
        ppp = repr(pp)
        if p == "dtype":
          try:
            ppp = pp.__name__
          except AttributeError:  
            ppp = str(pp)
        elif p == "flags":
          ppp = "\n" + ppp + "\n"
          if not pp.c_contiguous:
            pass
          elif not hasattr(pp, "owndata") or not pp.owndata:
            pass
          elif not hasattr(pp, "writeable") or pp.writeable:
            pass
          elif not hasattr(pp, "aligned") or pp.aligned:
            pass
          elif hasattr(pp, "updateifcopy") or pp.updateifcopy:  
            pass
          else:
            ppp = None
      if ppp is not None:   
        if not first: ret += ", "
        ret += p + " = " + ppp
        first = False
    ret += ")"
    return ret
  def get_data(self, view=True, refe=None):
    if refe is None: refe = self
    if self._parent is not None and self._index != 0:
      return self._parent.get_data(view, refe)
    data = self._data
    if data is None or view==False:
      return data
    assert data.flags.c_contiguous
    shape, dshape = refe.shape, data.shape
    if refe.dtype != data.dtype or dshape != shape:    
      if self._gpu:
        def dummy_alloc(nbytes):
          if nbytes > data.nbytes: 
            raise ViewError(nbytes, data.nbytes, shape, data.shape, refe.dtype, data.dtype)
          return data.ptr
        v = pycuda.gpuarray.GPUArray(shape, refe.dtype, dummy_alloc)        
        data = v
      else:
        try:
          v = np.ndarray(refe.shape, refe.dtype, buffer=data)
        except TypeError:
          raise ViewError(data.size, data.dtype, data.shape, refe.dtype, refe.shape)
        data = v
    return data
  def set_data(self, data):
    if self._parent is not None and self._index != 0: 
      return self._parent.set_data(data)
    if self.gpu:
      assert data is None or (has_pycuda and isinstance(data, pycuda.gpuarray.GPUArray))
      #TODO: check gpucontext??
    else:
      assert data is None or isinstance(data, np.ndarray)
    self._data = data
  def get_hash(self):
    return self._hash
  def set_hash(self, hash):
    self._hash = hash
  def is_assigned(self):
    return self._assigned or self._index > 0      
  def cache(self):        
    if not self._assigned: 
      assert self._parent._cache_active
      self._parent._cache_value = True
    assert self._cache_active  
    self._current()._cache_value = True
  def __getattr__(self, attr):
    if attr not in self._attrs:
      attr2 = "_" + attr
      if attr2 in self._attrs:
        attr = attr2
        return getattr(self, attr)
      else:
        raise AttributeError(attr)
    data = self.get_data(view=False)
    if data is None:
      raise AttributeError(attr)
    return getattr(data, attr)
  def __setattr__(self, attr, value):
    attr2 = "_" + attr    
    if attr2 in self._attrs:
      raise AttributeError("%s: Cannot re-define read-only attribute '%s'" % (self._tvname, attr))  
    if attr not in self._attrs:
      object.__setattr__(self, attr, value)
      return      
    data = self.get_data(view=False)
    if data is None or not hasattr(data, attr):
      if attr == "dtype":
        value = np.dtype(value)
      object.__setattr__(self, attr, value)
      return
    oldvalue = getattr(data, attr)

    if oldvalue != value:
      if attr == "shape":
        object.__setattr__(self, attr, value)
        return
      elif attr == "flags":
        #ignore for now
        return
      else:
        raise AttributeError("%s: Cannot re-define '%s' from '%s' to '%s'" % (self._tvname, attr, oldvalue, value))  
    object.__setattr__(self, attr, value)
    
  def _current(self):
    mode = get_mode()
    if mode == "direct":
      return self
    elif mode == "unroll":      
      if self._parent is None or not self._parent._assigned:
        return self
      else:
        return self._parent._children[-1]
  def _previous(self):
    assert self._index > 0
    return self._parent._children[self._index-1]    
  def _assign(self, accumulate):
    assert not self._sharded
    mode = get_mode()
    if mode == "direct":
      return self
    elif mode == "unroll":
      data = self.get_data()
      if data is not None:
        raise ValueError("TVArray '%s' was constructed with template data, assignment not (yet) supported" % self.name())
      if not self._assigned:
        assert not accumulate
        self._assigned = True
        return self
      else:
        #We must be the original (first) array in the series
        assert self._index == 0, self._index        
        return self._parent._assign(accumulate)
    else:
      raise ValueError(mode)
  def _shard(self, shapes, offsets, write):
    assert not self._sharded
    assert self.dtype is not None
    assert self.shape is not None
    if self.flags is not None:
      assert self.flags.c_contiguous
    if write:
      assert not self._join_assigned or self._index == 0
      if self._parent is None or not self._parent._assigned:
        pindex = None
      else:
        pindex = str(int(self._index)+1)
      self._sharded = "write"
      shardclass = TVWriteShard      
      
    else:
      self._sharded = "read"
      pindex = None
      assert self._parent._assigned or self._data is not None
      if self._assigned and self._index > 0:
        prevchild = self._parent._children[self._index-1]
        if prevchild._join_assigned:
          raise ValueError("Cannot define read shards for '%s', because write shards were defined for the previous assignment" % 
                           namestr(self.name()))
      shardclass = TVReadShard
    self._shards = []
    for n in range(len(offsets)):
      shard = shardclass(self, n, shapes[n], offsets[n], pindex=pindex)
      self._shards.append(shard)
  def join(self):   
    current = self._current()
    assert current._sharded
    assign = False
    if current._sharded == "write":
      for shard in current._shards:
        assert shard._assigned, namestr(shard.name())
      assign = True      
    current._sharded = False
    current._shards = None
    if assign:
      if not self._parent._assigned:      
        self._assign(accumulate=False)
      else:
        self._parent._assign(accumulate=False)
      self._join_assigned = True
    return self
  def rsplit(self):
    return self._split(False)  
  def wsplit(self):
    return self._split(True)  
  def _split(self, write):
    current = self._current()
    assert current.dtype is not None
    assert current.shape is not None
    assert len(current.shape) > 1
    strides = current.get_strides()
    shapes = [current.shape[1:] for n in range(current.shape[0])]
    offsets = [n*strides[0] for n in range(current.shape[0])]    
    current._shard(shapes, offsets, write)
    ret = []
    ret[:] = current._shards
    return ret
  def rchunks(self, chunksize):
    return self._chunks(chunksize, False)
  def wchunks(self, chunksize):
    return self._chunks(chunksize, True)  
  def _chunks(self, chunksize, write):
    current = self._current()
    assert current.dtype is not None
    assert current.shape is not None
    strides = current.get_strides()
    delims = list(range(0, current.shape[0], chunksize))
    if delims[-1] != current.shape[0]:
      delims.append(current.shape[0])
    shapes = [(delims[n+1]-delims[n],)+current.shape[1:] for n in range(len(delims)-1)]
    offsets = [n*strides[0] for n in range(0, current.shape[0], chunksize)]        
    current._shard(shapes, offsets, write)
    ret = []
    ret[:] = current._shards
    return ret
  def __len__(self):
    return self.shape[0]  
  def dict(self, state):
    d = {}
    for attr in self._attrs + ("_index", "_tvname", "_cache_value", "_cache_active", "_assigned"):
      v = getattr(self, attr)
      d[attr] = v
    if state is not None and self._parent is not None:
      parentindex = state.find_data_container(self._parent)
      d["_parent"] = parentindex
    return d
  @classmethod
  def fromdict(cls, dic, state):
    self = object.__new__(cls)    
    self._cache = TVCache(self) 
    for attr in self._attrs + ("_index", "_tvname", "_cache_value", "_cache_active", "_assigned"):
      if attr in dic:
        setattr(self, attr, dic[attr])
    if "_parent" in dic:
      self._parent = state.data_containers[dic["_parent"]]
      assert self._parent is not None, dic["_parent"]
      for n in range(len(self._parent._children), self._index+1):
        self._parent._children.append(None)
      self._parent._children[self._index] = self
    self.set_hash(None)
    return self

class TVBaseShard(TVData):
  _accessor = None
  _cache = None
  def __init__(self, parent, index, shape, offset, pindex=None):
    
    self._init = True
    self._parent = parent
    self._index = index
    self._offset = offset
    self.shape = shape
    self.flags = parent.flags
    self.dtype = np.dtype(parent.dtype)
    self._gpu = parent._gpu
    self._gpucontext = parent._gpucontext
    self._cache_active = False
    self.check_bounds()
    self._init = False
    
    self._cache = TVCache(self) 
    if get_mode() == "unroll":
      register(self)   
  
  def __setattr__(self, attr, v):
    if attr == "_init":
      return object.__setattr__(self, attr, v)
    if not self._init and attr in ("shape", "dtype"):
      if attr == "dtype":
        v = np.dtype(v)
      vv = getattr(self, attr)
      if v != vv:
        raise ValueError("Cannot change shard read-only property '%s' from '%s' to '%s'" % (attr, vv, v))
    return TVData.__setattr__(self, attr, v)  

  def check_bounds(self):  
    itemsize = np.dtype(self.dtype).itemsize
    nb = itemsize
    for s in self.shape: nb *= s
    nb_parent = itemsize
    for s in self._parent.shape: nb_parent *= s
    assert nb + self._offset <= nb_parent, (self._parent.name(), self._index, nb, self._offset, nb_parent, self.dtype, itemsize, self.shape, self._parent.shape)
  
  def get_data(self, view=True):    
    if get_mode() == "direct":
      assert self._parent._sharded
    data = self._parent.get_data(view=False)
    if view == False:
      return data
    assert data is not None
    self.check_bounds()
    if self._gpu:
      assert has_pycuda
      ptr = np.intp(data.ptr) + self._offset #raw pointer!
      d = pycuda.gpuarray.GPUArray(self.shape, self.dtype, lambda nbytes: ptr)
    else:
      d = np.ndarray(self.shape, buffer=data, offset=self._offset, dtype=self.dtype)
    return d  
  def _current(self):
    return self
  def params(self, hashing):
    return self._parent.params(hashing)
  def __getattr__(self, attr):
    if attr in ("gpu", "gpucontext"):
      return getattr(self, "_"+attr)
    else:
      raise AttributeError(attr)
  def cache(self):        
    raise TypeError("Shards cannot be cached")
  def __len__(self):
    return self.shape[0]
  def __eq__(self, other):
    return self._parent == other
  def __ne__(self, other):
    return not self.__eq__(other)

class TVWriteShard(TVBaseShard):
  def __init__(self, *args, **kwargs):
    TVBaseShard.__init__(self, *args, **kwargs)
    self._assigned = False    
    self._pindex = kwargs["pindex"]
  def name(self):
    p = self._pindex
    if p is None and len(self._parent._parent._children) > 1:
      p = '0'
    return (self._parent.name()[0], p, str(self._index))    
  def _assign(self, accumulate):
    assert not accumulate
    assert not self._assigned
    assert self._parent._sharded == "write"
    self._assigned = True
    return self
  def dict(self, state):
    d = {}
    attrs = ("assigned", "index", "offset")
    for attr in attrs:
      d[attr] = getattr(self, "_" + attr)      
    d["shape"] = self.shape
    d["parent"] = state.find_data(self._parent)
    d["pindex"] = self._pindex
    d["write"] = True
    return d
  @classmethod
  def fromdict(cls, dic, state):
    parent = state.data[dic["parent"]]
    index, shape, offset = dic["index"], dic["shape"], dic["offset"]
    self = cls(parent, index, shape, offset, pindex=dic["pindex"])
    self._assigned = dic["assigned"]
    return self

class TVReadShard(TVBaseShard):
  _offset = None
  def name(self):
    return self._parent.name() + (str(self._index),)
  def _assign(self, accumulate):
    raise TypeError("Read shards cannot be assigned to")
  def dict(self, state):
    d = {}
    attrs = ("index", "offset")
    for attr in attrs:
      d[attr] = getattr(self, "_" + attr)      
    d["shape"] = self.shape
    d["parent"] = state.find_data(self._parent)
    d["write"] = False
    return d
  @classmethod
  def fromdict(cls, dic, state):
    parent = state.data[dic["parent"]]
    index, shape, offset = dic["index"], dic["shape"], dic["offset"]
    self = cls(parent, index, shape, offset)
    return self
  def get_hash(self):
    return self._hash
  def search_cache_hash(self):    
    from .TVCache import TVHash
    p = self._parent
    assert p._cache_active  
    if self._parent.get_data() is not None and self._parent._cache_active and not self._parent._cache._readshard_saved:      
      d = self._parent.get_data()      
      arrayhash = self._parent.get_hash()
      if arrayhash is None:
        arrayhash = TVHash(d)        
        self._parent.set_hash(arrayhash)
      self._parent._cache.arrayhash = arrayhash
      self._parent._cache.save_readshard_cache(arrayhash)    
      h = p._cache.search_readshard(self._offset)
      assert h is not None, self.name()
    else:
      h = p._cache.search_readshard(self._offset)
    self._hash = h
  
  
from .TVCache import TVCache
from .TVContext import register
