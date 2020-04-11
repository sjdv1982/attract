from .Space import Space, _logger as logger
from . import Handler
from .Handler import UniNullHandler, UniMethod, HandlerError, Constructor
from weakref import WeakValueDictionary

def decompose_attribute(attr):
  ret = []
  curr = attr
  try:
    while 1:
      pdot, phook = curr.find("."), curr.find("[")
      if pdot > -1 and (phook == -1 or pdot < phook):
        if pdot > 0:
          ret.append(curr[:pdot])
        curr = curr[pdot+1:]
      elif phook > -1:
        phook2 = curr.index("]")
        ret.append(curr[:phook])
        index = int(curr[phook+1:phook2])
        ret.append(index)
        curr = curr[phook2+1:]
      else:
        break
  except:
    raise ValueError("Invalid attribute expression", attr)
  if len(curr): ret.append(curr)
  return tuple(ret)

class BaseNode(object):
  _handlers = []
  @classmethod
  def addHandler(cls, handler):
    assert isinstance(handler, Handler.Handler)
    cls._handlers.insert(0, handler) #LIFO
  @classmethod    
  def removeHandler(cls, handler):
    assert isinstance(handler, Handler)
    cls._handlers.remove(handler)
  @classmethod    
  def getHandlers(cls):
    return cls._handlers
  
  def __init__(self, nodetype, root, parent, name):
    self.nodetype = nodetype
    if root is None and parent is None: root = self
    self.root = root
    self.parent = parent
    self.name = name
    self.space = Space()
    self._children = {}
    self._links = WeakValueDictionary()
    self._children_names = []
    self._myhandlers = None
    self._lock_children = False
    self.currentAction = None
    self.assigned = list(self.__dict__.keys())
  
  def __setattr__(self, attr, value):
    if hasattr(self, "assigned"):
      if attr not in self.assigned: raise AttributeError(attr)
    self.__dict__[attr] = value
    
  def getChildren(self):
    return self._children_names

  def getLinks(self):
    return self._links.keys()
    
  def hasChild(self, name):
    return name in self._children
    
  def addChild(self, nodetype, name):
    self.currentAction = "addChild"
    logger("addChild", self.nodetype, self.name, nodetype, name)
    assert not self._lock_children
    assert name not in self._children, name
    self._lock_children = True
    child = None
    try:
      child = Node(nodetype, self.root, self, name)
      child.sendMessage("init")
      self._children[name] = child
      self._children_names.append(name)      
    finally:
      self._lock_children = False  
    child.sendMessage("post-init")  
    self.currentAction = None
    return child

  def attachChild(self, child, name, position=None):
    self.currentAction = "attachChild"
    logger("attachChild", self.nodetype, self.name, name, child.nodetype, position)
    assert not self._lock_children
    assert name not in self._children, name
    self._lock_children = True
    try:
      self._children[name] = child
      if position is None:
        self._children_names.append(name)      
      else:
        self._children_names.insert(position, name)
        for n in range(position+1, len(self._children)):
          childname = self._children_names[n]
          self._children[childname].sendMessage("child-index-increase")                    
    finally:
      self._lock_children = False  
    child.parent = self
    child.name = name
    set_root(child, self.root)
    child.sendMessage("attach")
    self.currentAction = None
    return child

  def attachLink(self, link, name):
    self.currentAction = "attachLink"
    logger("attachLink", self.nodetype, self.name, name, link.nodetype)
    assert name not in self._children, name
    assert name not in self._links, name
    self._links[name] = link
    self.currentAction = None
    return link
    
  def insertChild(self, nodetype, name, position):
    self.currentAction = "insertChild"
    logger("insertChild", self.nodetype, self.name, nodetype, name, position)
    assert not self._lock_children
    assert name not in self._children, name
    self._lock_children = True
    child = None
    try:
      child = Node(nodetype, self.root, self, name)
      child.sendMessage("init")
      self._children[name] = child
      self._children_names.insert(position, name)
      for n in range(position+1, len(self._children)):
        childname = self._children_names[n]
        self._children[childname].sendMessage("child-index-increase")            
    finally:
      self._lock_children = False  
    child.sendMessage("post-init")    
    self.currentAction = None
    return child
    
  def removeChild(self, name):
    self.currentAction = "removeChild"
    logger("removeChild", self.nodetype, self.name, name)
    assert not self._lock_children
    assert name in self._children, name
    self._lock_children = True
    child = None
    try:    
      position = self._children_names.index(name)
      child = self._children[name]
      child.sendMessage("destroy")
      self._children_names.remove(name)      
      self._children.pop(name)
      for n in range(position, len(self._children)):
        childname = self._children_names[n]
        self._children[childname].sendMessage("child-index-decrease")
    finally:
      self._lock_children = False 
    child.sendMessage("post-destroy") 
    self.currentAction = None

  def detachChild(self, name):  
    self.currentAction = "detachChild"
    logger("detachChild", self.nodetype, self.name, name)
    assert not self._lock_children
    assert name in self._children, name
    self._lock_children = True
    child = None
    try:    
      position = self._children_names.index(name)
      child = self._children[name]
      child.sendMessage("pre-detach")
      child.parent = None
      child.name = None
      self._children_names.remove(name)      
      self._children.pop(name)
      for n in range(position, len(self._children)):
        childname = self._children_names[n]
        self._children[childname].sendMessage("child-index-decrease")
      child.sendMessage("post-detach")
    finally:
      self._lock_children = False       
      self.currentAction = None
    return child

  def detachAndReplaceChild(self, name, newchild):  
    self.currentAction = "detachAndReplaceChild"
    logger("detachAndReplaceChild", self.nodetype, self.name, name, newchild.nodetype)
    assert not self._lock_children
    assert newchild.parent is None
    assert newchild.name is None
    assert name in self._children, name
    self._lock_children = True
    child = None
    try:    
      position = self._children_names.index(name)
      child = self._children[name]
      child.sendMessage("pre-detach")
      child.parent = None
      child.name = None
      self._children[name] = newchild
      newchild.parent = self
      newchild.name = name
      set_root(newchild, self.root)
      newchild.sendMessage("attach")
      child.sendMessage("post-detach")
    finally:
      self._lock_children = False 
      self.currentAction = None    
    return child

  def replaceChild(self, name, newchild):  
    self.currentAction = "replaceChild"
    logger("replaceChild", self.nodetype, self.name, name, newchild.nodetype)
    assert not self._lock_children
    assert newchild.parent is None
    assert newchild.name is None
    assert name in self._children, name
    self._lock_children = True
    child = None
    try:    
      position = self._children_names.index(name)
      child = self._children[name]
      child.sendMessage("destroy")
      self._children[name] = newchild
      newchild.parent = self
      newchild.name = name
      set_root(newchild, self.root)
      newchild.sendMessage("attach")
      child.sendMessage("post-destroy")
    finally:
      self._lock_children = False 
      self.currentAction = None
    return child
    
  def sendMessage(self, message, args = ()):
    self.currentAction = "sendMessage"
    logger("sendMessage", self.nodetype, self.name, message, args)
    if self._myhandlers is None:
      self._buildMyHandlers()
    try:
      matched = False
      while 1:
        matches = []
        fallbacks = []
        
        claimed = False
        for handler in self._myhandlers:
          if handler.nodetype is not None and handler.nodetype != self.nodetype: continue #TODO: see buildMyHandlers
          if handler.message is not None and handler.message != message: continue
          poll = handler.poll(message, self.nodetype)
          if poll == Handler.CLAIM:
            result = handler.invoke(message, self, args)
            assert result != Handler.NO_MATCH
            claimed = True
            break
          elif poll == Handler.MATCH:    
            matches.append(handler)
          elif poll == Handler.NO_MATCH:  
            continue          
          elif poll == Handler.FALLBACK:    
            fallbacks.append(handler)
          else:
            raise Exception("Unknown handler poll result", handler, result)
        
        if claimed: 
          matched = True
          break
                
        for handler in matches:
          result = handler.invoke(message, self, args)
          if result == Handler.CLAIM: 
            matched = True
            break
          elif result == Handler.MATCH:   
            matched = True
            continue
          elif result == Handler.NO_MATCH: 
            continue
          else:
            raise HandlerError("Unknown handler invoke result", handler._invokefunc, result)
        
        if matched: break

        for handler in fallbacks:
          result = handler.invoke(message, self, args)
          if result == Handler.CLAIM: 
            matched = True
            break
          elif result == Handler.MATCH:   
            matched = True
            continue
          elif result == Handler.NO_MATCH: 
            continue
          else:
            raise HandlerError("Unknown handler invoke result", handler, result)
        
        break
        
      if not matched:
        raise HandlerError("No handler for this message", message, self.nodetype, self.name )
    finally:
      self.currentAction = None
      
  def morph(self, newnodetype):
    self.currentAction = "morph"
    logger("morph", self.nodetype, self.name, newnodetype)
    self._myhandlers = None
    oldnodetype = self.nodetype
    self.nodetype = newnodetype
    try:
      for childname in list(self._children.keys()):
        if childname not in self._children: continue #may have been removed in the meantime...
        self._children[childname].sendMessage("parent-morph", (oldnodetype, newnodetype))
      if self.parent is not None:
        self.parent.sendMessage("child-morph", (self.name, oldnodetype, newnodetype))
    finally:
      self.currentAction = None
      
  def _buildMyHandlers(self):
    myhandlers = []
    nodetype = self.nodetype
    """
    for handler in self._handlers:
      if handler.nodetype is not None and handler.nodetype != nodetype: continue
      myhandlers.append(handler)
    self._myhandlers = myhandlers  
    """
    self._myhandlers = self._handlers

  def __getitem__(self, item):
    try:
      return self._children[item]
    except KeyError:
      if item not in self._links: raise
      return self._links[item]
      
    
class Node(BaseNode): pass
Node.addHandler(UniNullHandler("init"))
Node.addHandler(UniNullHandler("parent-morph"))
Node.addHandler(UniNullHandler("child-morph"))
Node.addHandler(UniNullHandler("child-index-increase"))
Node.addHandler(UniNullHandler("child-index-decrease"))
Node.addHandler(UniNullHandler("post-init"))
Node.addHandler(UniNullHandler("attach"))
Node.addHandler(UniNullHandler("destroy"))
Node.addHandler(UniNullHandler("post-destroy"))
Node.addHandler(UniNullHandler("pre-detach"))
Node.addHandler(UniNullHandler("post-detach"))

def reparent(message, node, attribute, newparent, newname, position = None, replace = False, keep_link = False):
  assert isinstance(attribute, tuple), attribute
  assert len(attribute)
  for a in attribute:
    assert isinstance(a, str) or isinstance(a, int), a
  childname = attribute[0]
  if len(attribute) == 1:    
    childname = attribute[0]
    child = node.detachChild(childname)
    if keep_link == True:
      node.attachLink(child, childname)
    if replace == True:
      newparent.attachChild(child, newname, position)
    else:  
      newparent.replaceChild(newname, child)
  else:
    assert childname in node.getChildren() or childname in node.getLinks(), (node.getChildren(), node.getLinks())
    node[childname].sendMessage(message, args = (attribute[1:], newparent,newname, position, replace, keep_link) )

def set_root(node, root):
  node.root = root
  for child in node.getChildren():
    set_root(node[child], root)
    
import functools    
Node.addHandler(UniMethod("reparent", functools.partial(reparent, "reparent")))
Node.addHandler( Constructor("dummy", lambda node: None) )

class Root(Node):
  def __init__(self):
    self.root = self
    Node.__init__(self, "root", None, None, None)
    
    
