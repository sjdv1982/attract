unimethodlist = []
unimethods = {}
macros = set()
statements = {}
converters = []
inheritedtypes = {}

class general(object):
  pass

def defineunimethod(name, pointer):
  from . import corelock
  if corelock:
    raise Exception("system error: cannot define new unimethod, .spy files have already been compiled")
  if name in unimethods:
    raise Exception("system error: duplicate definition of unimethod %s" % name)
  unimethods[name] = pointer
  unimethodlist.append(name)

def definemacro(pointer):
  from . import corelock
  if corelock:
    raise Exception("system error: cannot define new macro, .spy files have already been compiled")
  macros.add(pointer)

def definestatement(name, pointer):
  from . import corelock
  if corelock:
    raise Exception("system error: cannot define new statement, .spy files have already been compiled")
  if name in statements:
    raise Exception("system error: duplicate definition of statement %s" % name)
  statements[name] = pointer  
