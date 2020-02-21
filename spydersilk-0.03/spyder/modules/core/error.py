error = {} #spyder.modules.core.error

import re, ast
quotematch = re.compile(r'(([\"\']).*?\2)')
curlymatch = re.compile(r'{[^{}]*?}')

def parse_error(s):
  ret = []
  pq = list(quotematch.finditer(s) )
  currpos = 0
  for pnr in range(0, len(pq), 2):    
    p1, p2 = pq[pnr], pq[pnr+1]
    mid = s[p1.end():p2.start()].strip().replace("\n","")
    s1, s2 = p1.group(0).strip(), p2.group(0).strip()
    if mid != "=>":
      raise ValueError("Malformed error statement: \n %s\n    %s\n %s\n'%s' should be '=>'" % (s1,mid, s2,mid))
    ss1 = ast.literal_eval(s1)
    ss2 = ast.literal_eval(s2)
    ret.append((ss1, ss2))
    currpos = p2.end()
  return ret
  
def _update_error(cls, statement, message):  
  t = cls.typename()
  if t not in error: return
  e = error[t]
  if statement not in e: return
  if message is not None:
    for p in curlymatch.finditer(message):        
      expr = p.group(0)[1:-1]  
      spyder.safe_eval(expr)  
  e[statement] = message
  
def _register_error(cls, statement, message):
  spyder.safe_eval(statement)
  if message is not None:
    for p in curlymatch.finditer(message):        
      expr = p.group(0)[1:-1]  
      spyder.safe_eval(expr)
  error[cls.typename()][statement] = message
    
def retrieve_message(statement,typ):
  t = typ.typename()
  errpath = spyder.errorpath
  message = None
  while 1:      
    if t not in error: break
    errors = error[t]
    if statement not in errors or errors[statement] in ("", None): break
    message = errors[statement]
    break
  if message is None:
    return statement
  return message 
  
def _assert(self, expression, statement):
  import inspect, traceback
  lastframe = inspect.currentframe().f_back
  result = eval(expression, lastframe.f_globals, lastframe.f_locals)
  if isinstance(result, tuple) and len(result) == 2 and (result[0] is True or result[0] is False):
    if not result[0]:
      message = retrieve_message(statement, self)
      if message == statement: message += " => " + str(result[1])      
    else:
      return
  else:
    if not result:
      message = retrieve_message(statement, self)
    else:
      return
  message_expr = ""
  curpos = 0
  errclass = AssertionError
  if message != statement: errclass = ValidationError
  for p in curlymatch.finditer(message):    
    message_expr += message[curpos:p.start()]
    expr = p.group(0)[1:-1]
    try:
      expr_result = str(eval(expr, lastframe.f_globals, lastframe.f_locals))
    except:
      expr_result = "{" + "".join(traceback.format_exception_only(*sys.exc_info()[:2])).replace("\n","") + "}"
    message_expr += expr_result      
    curpos = p.end()
  message_expr += message[curpos:]
  raise errclass(message_expr)
    
    
def _raise(self, expression, statement):  
  #print("RAISE: %s, %s" % (self.typename(), expression))
  import inspect  
  lastframe = inspect.currentframe().f_back  
  exec(statement, lastframe.f_globals, lastframe.f_locals)
  #TODO
  #import sys; sys.exit()
  
def insert_error(line, err):
  assert err in ("raise", "assert"), err
  
  errorline = "spyder.core._register_error(cls, %s, None)" % repr(line.strip()) 
  pos = line.rindex(err)
  if err == "assert":
    line2 = line[pos+len("assert"):].strip()
    execline = line[:pos] + "spyder.core._assert(self, %s, %s)" % (repr(line2), repr(line.strip()))
  else:
    line2 = line[pos+len("raise"):].strip()
    execline = line[:pos] + "spyder.core._raise(self, %s, %s)" % (repr(line2), repr(line.strip()))
  return execline, errorline
  
def generate_error(typename, parentnames, source, members, deleted_members, block):  
  if block is None: return None, None
  errors = parse_error(block)
  if not len(errors): return None, None
  s = """  @classmethod
  def _update_errors(cls):
    pass
"""
  for message, statement in errors:
    s += "    spyder.core._update_error(cls, %s, %s)\n" % (repr(message),repr(statement))
  return s, None
  
spyder.defineunimethod("error", generate_error)  