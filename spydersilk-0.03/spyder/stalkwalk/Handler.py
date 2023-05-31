CLAIM = "claim"
MATCH = "match"
NO_MATCH = "no_match"
FALLBACK = "fallback"

class HandlerError(Exception): pass

class Handler:
  def __init__(self, pollfunc, invokefunc, nodetype = None, message = None):
    self.nodetype = nodetype
    self.message = message
    self._pollfunc = pollfunc
    self._invokefunc = invokefunc
  def poll(self, message, nodetype):
    return self._pollfunc(message, nodetype)
  def invoke(self, message, node, args):
    return self._invokefunc(message, node, args)    

def UniNullHandler(message):
  return Handler (
    pollfunc = lambda message, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: CLAIM,
    message = message
  )
    
def NullHandler(nodetype, message):
  return Handler (
    pollfunc = lambda message, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: CLAIM,
    nodetype = nodetype,
    message = message,
  )

def Constructor(nodetype, func):
  return Handler (
    pollfunc = lambda message, nodetype: CLAIM,
    invokefunc = lambda message, node, args: func(node),
    nodetype = nodetype,
    message = "init",
  )

def UniMethod(message, func):
  return Handler (
    pollfunc = lambda m, nodetype: CLAIM,
    invokefunc = lambda message, node, args: func(node, *args),    
    message = message,
  )
  
def Method(nodetype, message, func):
  return Handler (
    pollfunc = lambda m, nodetype: CLAIM,
    invokefunc = lambda message, node, args: func(node, *args),    
    nodetype = nodetype,
    message = message,
  )

def XMethod(nodetype, message, func):
  """
  Defines extended method, receiving the name of the message
  """  
  return Handler (
    pollfunc = lambda m, nodetype: CLAIM,
    invokefunc = lambda message, node, args: func(message, node, *args),    
    nodetype = nodetype,
    message = message,
  )

def CMethod(nodetype, message, func):
  """
  Defines conditional method
  """
  return Handler (
    pollfunc = lambda m, nodetype: MATCH,
    invokefunc = lambda message, node, args: func(node, *args),    
    nodetype = nodetype,
    message = message,
  )

def XCMethod(nodetype, message, func):
  """
  Defines conditional extended method, receiving the name of the message
  """  
  return Handler (
    pollfunc = lambda m, nodetype: MATCH,
    invokefunc = lambda message, node, args: func(message, node, *args),    
    nodetype = nodetype,
    message = message,
  )
  
def passdown(message, node, *args):
  for childname in list(node.getChildren()):
    if node.hasChild(childname): #could be removed now...
      child = node[childname]
      child.sendMessage(message, *args)
  return CLAIM
  
def PassDown(nodetype, message):
  return Handler (
    pollfunc = lambda message, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: passdown(message, node, *args),
    nodetype = nodetype,
    message = message,
  )

def passup(message, node, *args):
  if node.parent is not None:
    node.parent.sendMessage(message, *args)
  return CLAIM
  
def PassUp(nodetype, message):
  return Handler (
    pollfunc = lambda message, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: passup(message, node, *args),
    nodetype = nodetype,
    message = message,
  )

def UniPassUp(message):
  return Handler (
    pollfunc = lambda message, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: passup(message, node, *args),
    message = message,
  )
  
def FallBack(nodetype, message, func):
  return Handler (
    pollfunc = lambda m, nodetype: FALLBACK,
    invokefunc = lambda message, node, args: func(node, *args),    
    nodetype = nodetype,
    message = message,
  )
  