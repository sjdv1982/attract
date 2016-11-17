from . import namestr, TVArray, log
from .TVAllocator import Allocator
from .TVAction import TVTask, GPUAction
from .TVDepsGraph import TVDepsTask
from functools import partial
import sys
import numpy as np
import threading
import time
from collections import deque


class Evaluator(object):
  memallocator = None
  def evaluate_resource(self):
    raise NotImplementedError
  def run(self):
    raise NotImplementedError
  def init(self):
    pass

class SimpleEvaluator(Evaluator):    
  def __init__(self, graph):
    self.graph = graph
    self.done_resources = set()
    self.done_data = set()
    self.done_tasks = set()
    self.done_allocators = set()
    self.allocator = Allocator(self.graph, self)
    self.jobs = []

  def evaluate_allocator(self,allocator):
    if allocator in self.done_allocators: return
    d = self.graph.allocators[allocator]
    #print "EVALUATE ALLOCATOR", d
    self.jobs.append(partial(d.evaluate, self.allocator))
    self.done_allocators.add(allocator)

  def evaluate_data(self, data):
    if data in self.done_data: return
    #print "EVALUATE DATA", data
    #print
    d = self.graph.data[data]
    dd = d.tvdata()
    #print "DATA!", d    
    if d.resource:
      self.evaluate_resource(d.resource)
    if self.graph.cachedir is not None and dd._cache.cachekey is not None:
      arr = np.load(dd._cache.cachekey)
      dd.set_data(arr)      
    if d.task:
      self.evaluate_task(d.task)
      
    #print "/EVALUATE DATA", data
    #print
    self.done_data.add(data)
  
  def evaluate_task(self, task):
    def do_evaluate_task(task, tasknr, cache):
      log.log("EVALUATE", tasknr,  str(task).split("\n")[0].strip())
      #print "EVALUATE"
      #print task
      for c in task._cachekeys:
        log.log("LOAD " + namestr(c.name()), c._cache.cachekey)
      task._cachekeys_load()
      task.tvtask().evaluate(cache)
    if task in self.done_tasks: return
    #print "EVALUATE TASK", task
    #print
    d = self.graph.tasks[task]
    #print d
    for a in d.data:
      self.evaluate_data(a)
    for a in d.resources:
      self.evaluate_resource(a)
    if self.graph.cachedir is not None:
      self.jobs.append(partial(do_evaluate_task,d,task,cache=True))
    else:  
      self.jobs.append(partial(do_evaluate_task,d,task,cache=False))
    #print "/EVALUATE TASK", task
    #print
    self.done_tasks.add(task)
    
  def evaluate_resource(self, resource):
    if resource in self.done_resources: return
    #print "EVALUATE RESOURCE", resource
    #print
    d = self.graph.resources[resource]
    for a in d.data:
      self.evaluate_data(a)
    for a in d.resources:
      self.evaluate_resource(a)
    for a in d.allocators:
      self.evaluate_allocator(a)
      
    #print "/EVALUATE RESOURCE", resource
    #print
    self.done_resources.add(resource)

  def run(self):
    for j in self.jobs:
      j()


try:
  from threading import _Event as CPUEvent
except: #python3  
  from threading import Event as CPUEvent

try:
  import pycuda.driver
  from pycuda.autoinit import context
  class GPUEvent(object):
    _event = None
    def event(self):
      if self._event is None:
        self._event = pycuda.driver.Event(pycuda.driver.event_flags.BLOCKING_SYNC)
      return self._event      
except ImportError:
  class GPUEvent(object):
    def __init__(self):
      raise Exception("pycuda not found")

class GPUStream(object):
  """
  Delayed construction and execution of a CUDA stream
  Must be executed in the main thread!
  """
  def __init__(self, name, task, cache, gpu_events=[], gpu_out_events=[], cpu_out_events=[]):
    self.name = name
    self.task = task    
    self.events = gpu_events
    self.out_events = gpu_out_events
    self.cpu_out_events = cpu_out_events
    self.cache = cache
    self.running = False
  def runnable(self):
    return all([e._event is not None for e in self.events])
  def run(self):
    import pycuda.driver
    stream = pycuda.driver.Stream()             
    self.stream = stream
    self.waited_events = []
    for e in self.events:
      #print >> sys.stderr, "GPU EVENT WAIT", self.name, "for", e._eventname
      e = e.event()
      stream.wait_for_event(e)      
      ee = GPUEvent()
      self.waited_events.append(ee)
      ee.event().record(stream)    
    if self.task is not None: 
      log.log("EVALUATE GPU ", self.name, str(self.task).split("\n")[0].strip())
      self.task._cachekeys_load()
      tvtask = self.task.tvtask()
      try:
        tvtask.action.set_stream(stream) 
        tvtask.evaluate(self.cache)  #async
      finally:
        tvtask.action.set_stream(None)
    for e in self.out_events:
      e = e.event()      
      e.record(stream)    
    self.term_event = GPUEvent().event()
    self.term_event.record(stream)
    self.running = True
  def is_finished(self, start_event):
    if not self.stream.is_done: return False
    if not self.term_event.query(): return False
    for e,ee in zip(self.events, self.waited_events):
      eventname = e._eventname
      e = e.event()
      ee = ee.event()
      #print >> sys.stderr, "GPU EVENT WAITED", eventname, ee.time_since(start_event), e    
    #print >> sys.stderr, "GPU DONE", self.name, str(self.task).split("\n")[0].strip(), self.term_event.time_since(start_event)
    for e in self.out_events:
      eventname = e._eventname
      e = e.event()
      #print >> sys.stderr, "GPU EVENT SET", eventname, e.time_since(start_event)
    for e in self.cpu_out_events:
      #print >> sys.stderr, "CPU EVENT SET", e._eventname  
      e.set()
    return True
      

class CPUJob(object):
  def __init__(self, name, func, events=[], out_events=[]):
    self.name = name
    self.func = func
    self.events = events
    self.out_events = out_events
    for e in self.events:
      if not isinstance(e, (CPUEvent, GPUEvent)):
        raise TypeError(e)
    for e in self.out_events:
      if not isinstance(e, (CPUEvent, GPUEvent)):
        raise TypeError(e)      
  def run(self, parent):
    for e in self.events:
      if isinstance(e, CPUEvent):
        #print >> sys.stderr, "EVENT WAIT", e._eventname
        e.wait()
    
    gpu_in = [e for e in self.events if isinstance(e, GPUEvent)]
    gpu_out = [e for e in self.out_events if isinstance(e, GPUEvent)]
    if len(gpu_in) or len(gpu_out):
      stream_events = []
      if len(gpu_in):
        stream_event = CPUEvent()
        stream_event._eventname = "job gpu-in: " + self.name
        stream_events = [stream_event]
      gpustream = GPUStream(self.name, None, None, gpu_in, gpu_out, stream_events)
      parent.add_gpustream(gpustream)
    if len(gpu_in): 
      #print >> sys.stderr, "EVENT WAIT", stream_event._eventname
      self.stream_event = stream_event
      stream_event.wait()
            
    if self.func is not None:    
      self.func()
    for e in self.out_events:
      if isinstance(e, CPUEvent):
        #print >> sys.stderr, "EVENT SET", e._eventname
        e.set()
        
class CPUTask(object):
  def __init__(self, task, tasknr, cache, events=[], out_events=[]):
    self.task = task
    self.tasknr = tasknr
    self.name = str(task)
    assert isinstance(task, TVDepsTask), task
    assert not isinstance(task.tvtask().action, GPUAction)    
    self.cache = cache
    assert cache in (True, False), cache
    self.events = events
    self.out_events = out_events
    for e in self.events:
      if not isinstance(e, (CPUEvent, GPUEvent)):
        raise TypeError(e)
    for e in self.out_events:
      if not isinstance(e, CPUEvent):
        raise TypeError(e)      
  def run(self, parent):
    for e in self.events:
      if isinstance(e, CPUEvent):
        #print >> sys.stderr, "EVENT WAIT", e._eventname
        e.wait()

    gpu_in = [e for e in self.events if isinstance(e, GPUEvent)]
    if len(gpu_in):
      stream_event = CPUEvent()
      stream_event._eventname = "task gpu-in: %s" % str(self.task)
      gpustream = GPUStream(self.tasknr, None, None, gpu_in, [], [stream_event])
      parent.add_gpustream(gpustream)
      #print >> sys.stderr, "EVENT WAIT", stream_event._eventname
      stream_event.wait()
    
    log.log("EVALUATE CPU " + str(self.task).split("\n")[0].strip())
    self.task._cachekeys_load()
    tvtask = self.task.tvtask()
    tvtask.evaluate(self.cache)    

    for e in self.out_events:
      #print >> sys.stderr, "EVENT SET", e._eventname
      e.set()
  
class GPUTask(object):
  def __init__(self, task, tasknr, cache, events=[], out_events=[]):
    import pycuda
    self.task = task
    self.tasknr = tasknr
    self.name = str(task)
    assert isinstance(task, TVDepsTask), task    
    self.cache = cache
    assert cache in (True, False), cache    
    self.events = events
    self.out_events = out_events
    for e in self.events:
      if not isinstance(e, (CPUEvent, GPUEvent)):
        raise TypeError(e)
    for e in self.out_events:
      if not isinstance(e,(CPUEvent, GPUEvent)):
        raise TypeError(e)   
        
            
  def run(self, parent):
    for e in self.events:
      if isinstance(e, CPUEvent):
        #print >> sys.stderr, "EVENT WAIT", e._eventname
        e.wait()
        

    gpu_in = [e for e in self.events if isinstance(e, GPUEvent)]
    gpu_out = [e for e in self.out_events if isinstance(e, GPUEvent)]
    cpu_out = [e for e in self.out_events if isinstance(e, CPUEvent)]
    
    gpustream = GPUStream(self.tasknr, self.task, self.cache,  gpu_in, gpu_out, cpu_out)
    parent.add_gpustream(gpustream)
        
        
class ParallelEvaluator(Evaluator):   
  def __init__(self, graph):
    self.graph = graph
    self.done_resources = set()
    self.done_data = set()
    self.done_tasks = set()
    self.allocator = Allocator(self.graph, self)
    self.with_gpu = False
    self.jobs = []


  def evaluate_data(self, data):
    if data in self.done_data: return
    d = self.graph.data[data]
    dd = d.tvdata()
    events = []
    if d.resource:
      self.evaluate_resource(d.resource)
      event = self.events["resource", d.resource]
      events.append(event)      
    if self.graph.cachedir is not None and dd._cache.cachekey is not None:
      arr = np.load(dd._cache.cachekey)
      dd.set_data(arr)      
    if d.task:
      self.evaluate_task(d.task)
      event = self.events["task", d.task]
      events.append(event)            
    
    if d.task is not None and isinstance(self.graph.tasks[d.task].tvtask().action, GPUAction):
      self.with_gpu = True
      out_event = GPUEvent()
    elif d.tvdata()._gpu:
      self.with_gpu = True
      out_event = GPUEvent()
    else:  
      out_event = CPUEvent()
    out_event._eventname = "data: " + namestr(data)
    self.events["data", data] = out_event
    job = CPUJob("data job: " + namestr(data), None, events, [out_event])
    self.jobs.append(job)
    self.done_data.add(data)
  
  def evaluate_task(self, task):
    if task in self.done_tasks: return
    d = self.graph.tasks[task]
    events = []
    for a in d.data:
      self.evaluate_data(a)
      event = self.events["data", a]
      events.append(event)            
    for a in d.resources:
      self.evaluate_resource(a)
      event = self.events["resource", a]
      events.append(event)            
    if self.graph.cachedir is not None:
      cache = True
    else:  
      cache = False
    
    any_gpu = False    
    if isinstance(d.tvtask().action, GPUAction):
      any_gpu = True
    else:
      for a in d.tvtask()._outputs:
        if isinstance(a, TVArray) and a._gpu:
          any_gpu = True
    if not any_gpu:
      for a in d.data:
        data = self.graph.data[a]
        tvdata = data.tvdata() 
        if tvdata is not None and tvdata._gpu:
          any_gpu = True
          break
    if any_gpu:
      self.with_gpu = True
      out_event = GPUEvent()
      job = GPUTask(d, task, cache, events, [out_event])
    else:
      out_event = CPUEvent()            
      job = CPUTask(d, task, cache, events, [out_event])
    out_event._eventname = "task: " + str(task)
    self.events["task", task]  = out_event    
    self.jobs.append(job)
    self.done_tasks.add(task)
    
  def evaluate_resource(self, resource):
    if resource in self.done_resources: return
    d = self.graph.resources[resource]
    events = []
    any_cpu = False
    for a in d.data:
      self.evaluate_data(a)
      event = self.events["data", a]
      if isinstance(event, CPUEvent):
        any_cpu = True
      events.append(event)                  
    for a in d.resources:
      self.evaluate_resource(a)
      event = self.events["resource", a]
      if isinstance(event, CPUEvent):
        any_cpu = True      
      events.append(event)                  
      
    for a in d.allocators:
      allocator = self.graph.allocators[a]
      has_cache = False
      if allocator.alloc_command == "cache":
        assert not has_cache #Only one cache allocator
        f = functools.partial(allocator.evaluate, self.allocator)
        alloc_event = CPUEvent()
        job = CPUJob(a, f, events, [alloc_event])
        events = [alloc_event]
        any_cpu = True
        has_cache = True
      else:
        allocator.evaluate(self.allocator)

    if any_cpu or resource == "@return":
      out_event = CPUEvent()
      s = resource
      if resource != "@return":
        s = namestr(resource)
      out_event._eventname = "resource: " + s
      job = CPUJob("resource job: " + s, None, events, [out_event])
      self.jobs.append(job)      
    else:      
      out_event = GPUEvent()
      out_event._eventname = "gpu resource: " + namestr(resource)
      stream = GPUStream("gpu resource: " + namestr(resource), None, None, events, [out_event])
      self.add_gpustream(stream)
    self.events["resource", resource]  = out_event       
    self.done_resources.add(resource)

  def add_gpustream(self, stream):
    self._gpustreams.append(stream)
    #if threading.current_thread().name == 'MainThread':
    #  self.run_gpustreams()
    
  def run_gpustreams(self):
    if self.start_event is None:
      self.start_event = GPUEvent().event()
      self.start_event.record()    
      self.start_event.synchronize()
    
    while len(self._gpustreams):
      d = self._gpustreams.pop()          
      self._waiting_gpustreams.append(d)

    popped = []
    for dnr, d in enumerate(self._waiting_gpustreams):      
      if d.runnable(): 
        d.run()
        self._running_gpustreams.append(d)
        popped.append(dnr)
    for dnr in reversed(popped):
      self._waiting_gpustreams.pop(dnr)

    popped = []
    for dnr, d in enumerate(self._running_gpustreams):      
      if d.is_finished(self.start_event): 
        popped.append(dnr)
    for dnr in reversed(popped):
      self._running_gpustreams.pop(dnr)
    
    return (len(self._running_gpustreams) == 0)
    
  def run(self): 
    if self.with_gpu:
      import pycuda.tools
      self.memallocator = pycuda.tools.PageLockedMemoryPool()
    threads = deque()
    for j in self.jobs:
      t = threading.Thread(target=j.run, args=(self,))
      t.job = j
      threads.append(t)
      t.start()
    
    finished = False

    s_old = None    
    while not finished:
      tim = time.time()
      finished = self.run_gpustreams()
      if finished:
        #print "FINISHED?"
        while len(threads):
          t = threads.popleft()
          if t.is_alive():
            threads.appendleft(t)
            finished = False
            break        
      #sleep max. 0.1 ms  
      delay = 0.0001
      time.sleep(max(0, delay - (time.time() - tim)))
   
    self.memallocator.free_held()
    self.memallocator = None
    self.jobs = None
    self.events = None
   
  def init(self):
    self._gpustreams = deque()
    self.events = {}
    self._waiting_gpustreams = []
    self._running_gpustreams = []
    self.start_event = None
    
