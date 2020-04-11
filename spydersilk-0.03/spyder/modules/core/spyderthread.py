# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

import threading
from threading import Thread as tthread

class spythreadlist(threading.local):
  def __init__(self):
    self.l = []

class spythread(tthread):
  def __init__(self, group=None, target=None, name=None, args=(), kwargs={}):
    def threadfunc(*args,**kwargs):
      self.ret = target(*args, **kwargs)
    tthread.__init__(self,group,threadfunc,name, args,kwargs)
    self.args = args
    self.kwargs = kwargs
    self.setDaemon(True)
  def start(self):
    self.ret = None
    self.exception = None
    tthread.start(self)
  def run(self):
    try:
      tthread.run(self)
    except Exception as exc:
      self.exception = exc

threads = {}

def threadstart(target, args=(), kwargs={}, name=None):
  def joinname(name, nr):
    if nr == 0: return name
    else: return name + "-" + str(nr)
  if name == None:
    name = "Thread"
    nr = 1
  else:
    nr = 0
  while joinname(name, nr) in threads: nr += 1
  name = joinname(name, nr)
  threads[name] = spythread(target=target, name=name, args=args, kwargs=kwargs)
  threads[name].start()
  return name

def wait_all(*names):
  #This will wait for all of the threads to have finished
  #Any encountered errors will be raised immediately  
  import time
  ret = None
  try:
    while 1:
      finish_all = True
      for name in names:
        if name not in threads: raise Exception("Thread %s does not exist" % name)
        t = threads[name]
        if t == None: raise Exception("Thread %s has already returned" % name)
        e = t.exception      
        if e != None:
          raise type(e)("Exception in thread %s:\n%s" % (name, e.message ))
        if t.isAlive() == True: finish_all = False
      if finish_all == True:
        break
      time.sleep(1)
    ret = []
    for name in names: ret.append(threads[name].ret)    
  finally:
    for name in names: threads[name] = None
  if len(ret) == 1: return ret[0] 
  else: return tuple(ret)

def wait_any(*names):
  #This will wait for any of the threads to have finished
  #Any encountered errors will be raised immediately
  #Note that the remaining threads are NOT killed!
  import time
  ret = None
  while 1:
    finish_any = False
    for name in names:
      if name not in threads: raise Exception("Thread %s does not exist" % name)
      t = threads[name]
      if t == None: raise Exception("Thread %s has already returned" % name)        
      e = t.exception      
      if e != None:
        raise type(e)("Exception in thread %s:\n%s" % (name, e.message ))
      if t.isAlive() != True:
        finish_any = True
        break
    if finish_any == True:
      break
    time.sleep(1)
  threads[name] = None  
  return (name, t.ret)
        
