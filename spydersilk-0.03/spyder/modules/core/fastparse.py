# Copyright 2008-2014, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

def _fastparse_type(s, parsetype):
  dic = _fastparse(s)
  return parsetype.fromdict(dic)

def _fastparse(s):
  s = s.rstrip('\n').rstrip() + ","
  stack = [([],True),]
  curr,listmode = stack[-1]
  ident = 0
  objectlist = "ObjectList"
  for l in s.splitlines():
    if l.endswith("),"):
      if listmode is objectlist: #we are parsing an objectlist
        ll = l[ident+2*(curr[2]-1):]
        if ll == "),": #dedent
          if curr[2] > 0:
            curr[1] += l[ident:] + "\n"
            curr[2] -= 1
            if curr[2] == 0: #we are at outer level, parse what we have and reset
              txt = curr[1][:-2] #we have removed outer indentation from txt and we are effectively doing a 2nd pass here (suboptimal)
              curr[3].append(_fastparse_type(txt, spyder.__types__[curr[0]]))
              curr[:2] = None, None
          else: #dedent and leave objectlist mode
            curr[:] = curr[3]
            stack.pop()
            curr,listmode = stack[-1]
            ident -= 2
        else:
          if curr[2] == 0: #parsing an elemental value, e.g Float(1), at outer level
            ind = ll.index("(")
            typ = ll[2:ind]
            val = ll[ind+1:-2]
            curr[3].append(spyder.__types__[typ](val))
          else: #elemental value at inner level, treat it as continuation
            curr[1] += l[ident:] + "\n"
      else: #no objectlist, dedent
        stack.pop()
        curr,listmode = stack[-1]
        ident -= 2
    elif l.endswith(","): #continuation
      if listmode == True: #we're in an array, we expect unnamed items
        v = l[ident:-1]
        curr.append(v)      
      elif listmode is objectlist: #we're parsing an objectlist, gather it for later
        curr[1] += l[ident:] + "\n"       
      else: #we're in a class, we expect named items
        eq = l.index("=")
        k = l[ident:eq-1]
        v = l[eq+2:-1]
        curr[k] = v
    else: #endswith ( => indent
      if listmode == False: #we're in a class, we expect named items
        k = l[ident:l.find("=")-1]
        if l.find("ObjectList") > -1: #enter objectlist mode for the inner level
          new = [None, None, 0, []]
          listmode = objectlist
        elif l.endswith("Array ("): #enter array mode for the inner level
          new = []
          listmode = True
        else: #enter class mode for the inner level
          new = {}
          listmode = False          
        stack.append((new,listmode)) #push the stack
        curr[k] = new
      elif listmode is objectlist: #parsing objectlist is just gathering text
        if curr[2] == 0: #reset the text if we are at the outer level
          curr[0] = l[ident:-2]
          curr[1] = ""
        curr[1] += l[ident:] + "\n"
        curr[2] += 1
        ident -= 2 #to offset the +2 below, we don't want to increase ident read
      else:#we're in array mode,
        if l.find("ObjectList") > -1: #enter objectlist mode for the inner level
          new = [None, None, 0, []]
          listmode = "ObjectList"
        elif l.endswith("Array ("): #enter array mode for the inner level
          new = []
          listmode = True
        else: #enter class mode for the inner level
          new = {}
          listmode = False
        stack.append((new,listmode))
        curr.append(new)
      ident += 2 #increase indentation read
      curr = new
  assert len(stack) == 1 #when we're done, the stack must have been popped  
  return curr[0]

def fastparse(s):
  if bytes is not str and hasattr(s, "decode"): 
    s = s.decode()    
  try:
    parsetypename = s.splitlines()[0].split()[0]
  except Exception as e:
    e.__context__ = None
    raise spyder.core.ParseError("Object '%s' is unparsable" % type(s).__name__)    
  try:
    dic = _fastparse(s)
    return parsetypename, dic
  except Exception as e:
    e.__context__ = None
    raise spyder.core.ParseError("Object '%s' is unparsable" % type(s).__name__) 

