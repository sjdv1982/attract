import sys, os, numpy, struct, random, weakref, copy
import threading
from collections import deque, namedtuple
import atexit
import cPickle as pickle
import weakref

from . import namestr
from .hashing import Hasher, NumpyHasher

cachewritelock = threading.Lock()

class CacheIndex(object):
  def __init__(self, dir):   
    self.dir = dir
    self.indexfile = dir + os.sep + "index.pickle"
    self.data = {}
    self.datakeys = deque()
    self.indices = deque()
    self._written = False
    if os.path.exists(self.indexfile):
      self.data = pickle.load(open(self.indexfile))
      self.indices = deque([v[0] for v in self.data.values() if v[0] is not None])
      self.datakeys = deque([k for k in self.data.keys() if self.data[k][0] is not None])
      #print len(self.data), len(self.indices), len(self.datakeys)
  def search(self, outputindex, hashes):
    """
    print "SEARCH", self.dir, outputindex
    if (outputindex, hashes) not in self.data:
      print "MISS", self.dir, outputindex, hashes
      for k in self.data.keys():
        print "  DATA", k
    """
    return self.data.get((outputindex, hashes), None)
  def add(self, outputindex, hashes, array, cache_value):    
    if (outputindex, hashes) in self.data: 
      cachekey, arrayhash = self.data[(outputindex, hashes)]
      if cachekey is not None or not cache_value: 
        return cachekey, arrayhash
    assert isinstance(array, numpy.ndarray), type(array)        
    arrayhash = TVHash(array)
    if cache_value:
      cachekey = hash((outputindex, hashes))      
      cachekey = hex(cachekey)
      self.indices.append(cachekey)      
      self.data[outputindex, hashes] = cachekey, arrayhash
      self.datakeys.append((outputindex, hashes))
      self._written = False
      if not os.path.exists(self.dir):
        os.mkdir(self.dir)
      numpy.save(self.dir + os.sep + cachekey, array)   
    elif (outputindex, hashes) not in self.data:
      self.data[outputindex, hashes] = None, arrayhash      
      cachekey = None
      self._written = False  
    return cachekey, arrayhash
  def write_indexfile(self):
    #return ###
    with cachewritelock:
      newdata = self.data
      if os.path.exists(self.indexfile):
        newdata2 = pickle.load(open(self.indexfile))
        newdata2keys = set([k for k in newdata2.keys() if newdata2[k][0] is not None])
        if not newdata2keys <= set(self.datakeys):
          self._written = False
          newdata = newdata2
          newdata.update(self.data)
        write = ((newdata2keys != set(self.datakeys)) or set(newdata2.keys()) != set(newdata.keys())) and not self._written
      else:
        write = (len(newdata) and not self._written)
      if write:
        if not os.path.exists(self.dir):
          os.mkdir(self.dir)            
        f = open(self.indexfile, "w")
        pickle.dump(newdata, f)
        f.close()
      self._written = True  
  def __del__(self):    
    self.write_indexfile()

cache_indices = weakref.WeakValueDictionary()

def write_indexfiles():
  if cache_indices is None:
    return
  for k,v in cache_indices.items():
    v.write_indexfile()
    
atexit.register(write_indexfiles)

class TVHash(object):
   def __init__(self, array, hash=None):     
     self.shape = array.shape
     self.dtype = numpy.dtype(array.dtype).type
     self.hash = hash
     if hash is None:
       assert isinstance(array, numpy.ndarray)
       a = array
       if a.dtype == "float32":
         v = struct.unpack('i', '\x00\xFF\xFF\xFF')
         a = (a.view(dtype="int32") & v).view(dtype="float32")[:len(a)]       
       self.hash = NumpyHasher("sha1").hash(a)     
     else:
       raise Exception
   def __eq__(self, other):
     if not isinstance(other, TVHash): 
       return False
     return self.hash == other.hash and self.shape == other.shape and self.dtype == other.dtype
   def __ne__(self, other):
     return not self.__eq__(other)

   def __str__(self):
     return str((self.dtype, self.shape, self.hash))
   def __repr__(self):
     return str(self)
   def __hash__(self):
     return hash((self.dtype, self.shape, self.hash))

class TVCache(object):
  _readshard_saved = False
  def __init__(self, tvarray):    
    self.tvarray = weakref.ref(tvarray)
    self.graph = None
    self.searched = False
    self.cachekey = None
    self.arrayhash = None
        
  def search_join(self):
    assert self.graph is not None
    if self.searched: return (self.arrayhash is not None)
    outputindex, hashes = self._find_hashes_join()
    cacheindex = self.graph.cache_indices["::join"]
    result = cacheindex.search(outputindex, hashes)    
    self.searched = True    
    if result is None:      
      return False
    cachekey, arrayhash = result
    if cachekey is not None:
      cachekey = cacheindex.dir + os.sep + cachekey + ".npy"
    self.cachekey = cachekey
    self.arrayhash = arrayhash
    return True
  
  def search_readshard(self, offset):
    if self.arrayhash is None: 
      return None
    key = self.arrayhash, offset
    #print "SRCH", key
    indices = self.graph.cache_indices["::shard"]
    return indices.data.get(key, (None, None))[1]
  
  def search(self):
    assert self.graph is not None
    if self.searched: return (self.arrayhash is not None)
    task = self._find_task()
    #print "TASK", task
    result = self._find_hashes(task)
    outputindex, hashes = result
    cacheindex = self.graph.cache_indices[task.action.name()]
    result = cacheindex.search(outputindex, hashes)    
    self.searched = True    
    if result is None:      
      """
      print "CACHEMISS!", task
      from pprint import pprint
      pprint((outputindex, hashes))
      print
      for d in cacheindex.data:
        pprint(d)
      sys.exit()  
      """
      return False
    cachekey, arrayhash = result
    
    if cachekey is not None:
      #print >> sys.stderr, "HIT!", task, cachekey
      cachekey = cacheindex.dir + os.sep + cachekey + ".npy"
    else:
      #print >> sys.stderr, "PARTHIT!", task
      pass
    self.cachekey = cachekey
    self.arrayhash = arrayhash
    return True
 
  def _find_task(self):
    assert self.graph is not None    
    depsdata = self.graph.data[self.tvarray().name()]
    tasknr = depsdata.task
    task = self.graph.tasks[tasknr].tvtask()
    return task
  
  def _find_hashes(self, task):
    from .TVArray import TVBaseShard, TVReadShard, TVWriteShard
    outputindex = None
    for onr, o in enumerate(task._outputs):
      if o._cache is self:
        outputindex = onr
        break
    else:
      raise ValueError(self.tvarray().name(), str(task))    
    hashes = []
    for inr, i in enumerate(task._inputs):
      if isinstance(i, TVArray):
        readshard = None
        if isinstance(i, TVReadShard) and not i._parent._cache_active:
          readshard = i
          i = i._parent          
        if not i._cache_active and not isinstance(i, TVReadShard):
          d = self.graph.data[i.name()]
          assert d.hashtree is not None, namestr(i.name())
          hash = d.hashtree.bind(self.graph)
        elif i.get_hash() is not None:
          hash = i.get_hash()
        else:  
          assert not isinstance(i, TVReadShard), namestr(i.name())
          data = i.get_data()
          assert data is not None, (namestr(i.name()), str(task))
          hash = TVHash(data)
          i.set_hash(hash)
        if readshard is not None:
          hash = (readshard._offset, readshard.shape, hash)          
      else:
        hash = task.signature[inr]
      hashes.append(hash)
    hashes = tuple(hashes)
    return outputindex, hashes
  
  def _find_hashes_join(self):
    outputindex = -1
    hashes = []
    n = 0
    while 1:
      tvname2 = self.tvarray().name() + (n,)
      if tvname2 not in self.graph.data:
        break
      tasknr = self.graph.data[tvname2].task
      task = self.graph.tasks[tasknr]
      outputindex2, hashes2 = self._find_hashes(task)
      hashes.append((outputindex2, hashes2))
      n += 1
    hashes = tuple(hashes)
    return outputindex, hashes
    
  def save_join(self, array, cache_value):    
    assert self.graph.cache_indices is not None
    outputindex, hashes = self._find_hashes_join()

    taskname = "::join"    
    cachekey, arrayhash = self.graph.cache_indices[taskname].add(outputindex, hashes, array, cache_value)    
    return (outputindex, hashes), arrayhash

  def save(self, array, cache_value):    
    assert self.graph.cache_indices is not None
    task = self._find_task()
    outputindex, hashes = self._find_hashes(task)
    taskname = task.action.name()
    
    cachekey, arrayhash = self.graph.cache_indices[taskname].add(outputindex, hashes, array, cache_value)    
    #if cache_value:
    #  print  "SAVE!", self.tvarray().name(), cachekey, arrayhash
    self.save_readshard_cache(arrayhash)
    return (outputindex, hashes), arrayhash
  
  def save_readshard_cache(self, arrayhash):
    from .TVArray import TVReadShard
    a = self.tvarray()
    data = a.get_data()
    aname = a.name()
    assert data is not None and isinstance(data, numpy.ndarray), namestr(aname)
    self.readshard_cache = {}
    for dname in self.graph.data:
      if len(dname) < 3 or dname[:2] != aname: continue
      dd = self.graph.data[dname].tvdata()
      if not isinstance(dd, TVReadShard): continue
      ind = dname[2]
      ddata = dd.get_data()
      shardhash = TVHash(ddata) 
      key = (arrayhash, dd._offset)
      #print "SRC", dname[2], key
      self.graph.cache_indices["::shard"].data[key] = (None, shardhash)
    self._readshard_saved = True
    
from .TVArray import TVArray
      
    