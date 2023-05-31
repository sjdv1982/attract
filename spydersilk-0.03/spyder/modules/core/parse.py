# Copyright 2007-2013, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

bbytes = bytes
if bytes != str:
  import functools
  bbytes = functools.partial(bytes, encoding="utf8")

import re
parsetools = spyder.newnamespace()
parsetools.quotematch = re.compile(r'(?<![\\])(([\"\']).*?(?<![\\])\2)')
parsetools.brackmatch = re.compile(r'\[[^\[\]]*?\]')
parsetools.commamatch = re.compile(r'[ \n]*(,)[ \n]*(?=[^ \n]|$)')

class ParseError(TypeError, spyder.ConstructionError): pass

def listparse(s,mask, mode):
    if len(s) == 0: return []
    #mask all sub-lists
    mask0 = mask
    while 1:
        p = parsetools.brackmatch.finditer(mask0)
        pos = 0
        mask2 = ""
        for pp in p:            
            mask2 += mask0[pos:pp.start()] + "!" * (pp.end()-pp.start())
            pos = pp.end()
        mask2 += mask0[pos:]
        if pos == 0: break
        mask0 = mask2
    #split; we can't use split() because of the mask (we need start/end)
    elements = []
    p = parsetools.commamatch.finditer(mask2)
    pos = 0
    for pp in p:
        elements.append( [s[pos:pp.start()], mask[pos:pp.start()], mask2[pos:pp.start()]] )
        pos = pp.end(1)
    s = s.rstrip()
    if len(s[pos:].strip()) > 0: elements.append([s[pos:], mask[pos:len(s)], mask2[pos:len(s)]])

    values = []
    for enr, e in enumerate(elements):
      if e[0].startswith('\n'):
          e[0] = e[0][len('\n'):].lstrip()
          e[1] = e[1][len(e[1])-len(e[0]):]
      if len(e[0]) == 0:
        if enr < len(elements)-1: values.append("")
        continue
      v = e[0]
      m = e[1]
      v =v.lstrip()
      m = m[len(m)-len(v):]  
      bm = m.find('[')
      bm2 = m.rfind(']')      
      if v[0] in ('[','('): #list
          v = v[1:-1]
          m = m[1:-1]
          values.append(listparse(v,m,mode))
      elif v[0] == '{': #dict
          v = v[1:-1]
          m = m[1:-1]
          values.append(dictparse(v,m,mode))
      elif bm > -1: #function + args
          name = v[:bm].replace('\n','').strip()
          v = v[bm+1:bm2]
          m = m[bm+1:bm2]
          values.append(functionparse(name,v,m,mode))
      else: #value
          if (v[0] == "\"" and v[-1] == "\"") or (v[0] == "'" and v[-1] == "'"):
              values.append(v) #full quotes, for String and Data will deal with it
          elif v == "True": values.append(True)
          elif v == "False": values.append(False)                              
          else:
            try:
              if v.find(".") == -1: values.append(int(v))                
              else: values.append(float(v))
            except ValueError:
              values.append(v)
    return values

def dictparse(s,mask,mode):
    if len(s) == 0: return {}
    #mask all sub-lists
    mask0 = mask
    while 1:
        p = parsetools.brackmatch.finditer(mask0)
        pos = 0
        mask2 = ""
        for pp in p:
            mask2 += mask0[pos:pp.start()] + "!" * (pp.end()-pp.start())
            pos = pp.end()
        mask2 += mask0[pos:]
        if pos == 0: break
        mask0 = mask2
    #split; we can't use split() because of the mask (we need start/end)
    elements = []
    p = parsetools.commamatch.finditer(mask2)
    pos = 0
    for pp in p:
        elements.append( [s[pos:pp.start()], mask[pos:pp.start()], mask2[pos:pp.start()]] )
        pos = pp.end(1)
    s = s.rstrip()
    if len(s[pos:].strip()) > 0: elements.append([s[pos:], mask[pos:len(s)], mask2[pos:len(s)]])

    pairs = {}
    for e in enumerate(elements):
        if e[0].startswith('\n'):
            e[0] = e[0][len('\n'):].lstrip()
            e[1] = e[1][len(e[1])-len(e[0]):]
            e[2] = e[2][len(e[2])-len(e[0]):]
        if len(e[0]) == 0:
          continue
        v = e[0]
        m = e[1]
        m2 = e[2]
        v =v.lstrip()
        m = m[len(m)-len(v):]
        m2 = m2[len(m2)-len(v):]      
        bm = m.find('[')
        bm2 = m.rfind(']')
        ee = m2.find(":")
        if ee > -1: #pair
            k = v[:ee].rstrip()
            v = v[ee+1:]
            m = m[ee+1:]
            v =v.lstrip()
            m = m[len(m)-len(v):]
            if (k[0] != "\"" or k[-1] != "\"") and (k[0] != "'" or k[-1] != "'"):
                raise SyntaxError("Malformed dictionary key %s" % k)
            k = k[1:-1]
            if k in pairs:
                raise SyntaxError("Double key %s" % k)
            bm = m.find('[')
            bm2 = m.rfind(']')
            
            if v[0] in ('[','('): #list
                v = v[1:-1]
                m = m[1:-1]
                pairs[k] = listparse(v,m,mode)
            elif v[0] == '{': #dict
                v = v[1:-1]
                m = m[1:-1]
                pairs[k] = dictparse(v,m,mode)
            elif bm > -1: #function + args
                name = v[:bm].replace('\n','').strip()
                v = v[bm+1:bm2]
                m = m[bm+1:bm2]
                pairs[k] = functionparse(name,v,m,mode)
            else: #value
                if (v[0] == "\"" and v[-1] == "\"") or (v[0] == "'" and v[-1] == "'"):
                    pairs[k] = v #full quotes, for String and Data will deal with it
                elif v == "True": pairs[k] = True
                elif v == "False": pairs[k] = False
                else:
                  try:
                    if v.find(".") == -1: pairs[k] = int(v)                
                    else: pairs[k] = float(v)
                  except ValueError:
                    pairs[k] = v
        else:
            raise SyntaxError("Malformed dictionary pair %s" % e[0])
    return pairs

def functionparse(name,s,mask,mode):
    if len(s) == 0: return ([],{})    
    objectlistmode = False
    if name.find("ObjectList") > -1: objectlistmode = True
    #mask all sub-lists
    mask0 = mask
    while 1:
        p = parsetools.brackmatch.finditer(mask0)
        pos = 0
        mask2 = ""
        for pp in p:
            mask2 += mask0[pos:pp.start()] + "!" * (pp.end()-pp.start())
            pos = pp.end()
        mask2 += mask0[pos:]
        if pos == 0: break
        mask0 = mask2
    #split; we can't use split() because of the mask (we need start/end)
    elements = []
    p = parsetools.commamatch.finditer(mask2)
    pos = 0
    for pp in p:
        elements.append( [s[pos:pp.start()], mask[pos:pp.start()], mask2[pos:pp.start()]] )
        pos = pp.end(1)
    s = s.rstrip()
    if len(s[pos:].strip()) > 0: elements.append([s[pos:], mask[pos:len(s)], mask2[pos:len(s)]])
    values = []
    pairs = {}
    for e in elements:
        if e[0].startswith('\n'):
            e[0] = e[0][len('\n'):].lstrip()
            e[1] = e[1][len(e[1])-len(e[0]):]
            e[2] = e[2][len(e[2])-len(e[0]):]
        if len(e[0]) == 0:
          if enr < len(elements)-1: 
              if len(pairs) > 0:
                  raise SyntaxError("Unnamed argument after named one")                      
              values.append("")
          continue
        v = e[0]
        m = e[1]
        m2 = e[2]
        v =v.lstrip()
        m = m[len(m)-len(v):]
        m2 = m2[len(m2)-len(v):]      
        bm = m.find('[')
        bm2 = m.rfind(']')
        ee = m2.find("=")
        if ee > -1: #named value
            key = v[:ee].rstrip()
            v = v[ee+1:]
            m = e[1][ee+1:]
            v =v.lstrip()            
            m = m[len(m)-len(v):]
            bm = m.find('[')
            bm2 = m.rfind(']')                      
            if key in pairs:
                raise SyntaxError("Double key %s" % key)
            if v[0] in ('[','('): #list
                v = v[1:-1]
                m = m[1:-1]
                pairs[key] = listparse(v,m,mode)
            elif v[0] == '{': #dict
                v = v[1:-1]
                m = m[1:-1]
                pairs[key] = dictparse(v,m,mode)
            elif bm > -1: #function + args
                newname = v[:bm].replace('\n','').strip()
                v = v[bm+1:bm2]
                m = m[bm+1:bm2]
                pairs[key] = functionparse(newname,v,m,mode)
            else: #value
                if (v[0] == "\"" and v[-1] == "\"") or (v[0] == "'" and v[-1] == "'"):
                    pairs[key] = v #full quotes, for String and Data will deal with it
                elif v == "True": pairs[key] = True
                elif v == "False": pairs[key] = False
                else:
                  try:
                      if v.find(".") == -1: pairs[key] = int(v)                
                      else: pairs[key] = float(v)
                  except ValueError:
                      pairs[key] = v
                
        
        else:
            if len(pairs) > 0:
                raise SyntaxError("Unnamed argument after named one")            
            if v[0] in ('[','('): #list
                v = v[1:-1]
                m = m[1:-1]
                values.append(listparse(v,m,mode))
            elif v[0] == '{': #dict
                v = v[1:-1]
                m = m[1:-1]
                values.append(dictparse(v,m,mode))
            elif bm > -1: #function + args
                newname = v[:bm].replace('\n','').strip()
                v = v[bm+1:bm2]
                m = m[bm+1:bm2]
                o = functionparse(newname,v,m,mode)
                if objectlistmode and mode == "plain":
                    o = spyder.__types__[newname].fromdict(o)
                values.append(o)
            else:
                if (v[0] == "\"" and v[-1] == "\"") or (v[0] == "'" and v[-1] == "'"):
                    values.append(v) #full quotes, for String and Data will deal with it
                elif v == "True": values.append(True)
                elif v == "False": values.append(False)                
                else:
                  try:
                    if v.find(".") == -1: values.append(int(v))                
                    else: values.append(float(v))
                  except ValueError:
                    values.append(v)
    if mode == "Spyder":
      return spyder.__types__[name](*values,**pairs)      
    elif mode == "plain":
      if len(values) == 0: return pairs
      if len(pairs) == 0: return values
      return [values, pairs]

def parse(s, mode = "plain", return_tuple = False):
    if bytes is not str and hasattr(s, "decode"): 
      s = s.decode()    
    try:  
      p = parsetools.quotematch.finditer(s)
    except TypeError as e:
      e.__context__ = None
      raise spyder.core.ParseError("Object '%s' is unparsable" % type(s).__name__)
    pos = 0
    mask = ""
    for pp in p:
        mask += s[pos:pp.start()] + "*" * (pp.end()-pp.start())
        pos = pp.end()
    mask += s[pos:]
    mask = mask.replace('{', '[').replace('}',']')
    mask = mask.replace('(', '[').replace(')',']')#
    mask = mask.replace(r'\\', '**')
    v = s
    if len(v) == 0: return ((), {})
    m = mask
    bm = m.find('[')
    bm2 = m.rfind(']')
    cm = m.find(',')
    if v[0] in ('[','('): #list
        v = v[bm+1:bm2]
        m = m[bm+1:bm2]
        ret = listparse(v,m,mode)
        if return_tuple:
          return(ret, {})
        else: 
          return ret
    if v[0] == '{': #dict
        v = v[bm+1:bm2]
        m = m[bm+1:bm2]        
        ret = dictparse(v,m,mode)
        if return_tuple:
          return([], ret)
        else:
          return ret  
    elif bm > -1 and (cm == -1 or bm < cm): #function + args
        name = v[:bm].replace('\n','').strip()
        v = v[bm+1:bm2]
        m = m[bm+1:bm2]
        ret = functionparse(name,v,m,mode)
        if not return_tuple: return name, ret
        elif isinstance(ret, dict): return ([], ret)
        elif isinstance(ret, list): return (ret, {})
        else: return ret
    elif cm > -1:
        ret = listparse(v,m,mode)
        if not return_tuple: return ret
        return (ret,{})
    else: #value
        if not return_tuple: return v
        else: raise ParseError("Invalid value to parse: %s" % (str(s)))
       
_fromfilecode = """  @classmethod
  def fromfile(c, filename, fastparse=False):
    \"\"\"Auto-generated by Spyder:
     module core
     file parse.py
     function generate_parse   
    fromfile constructor\"\"\"
    n = spyder.linesep
    exc = ""
    try:
      try: #first attempt a parsing read, does not work with all data models
        s = spyder._file_load(filename).read()
        if fastparse:
          data = spyder.core.fastparse(s)
        else:
          data = spyder.core.parse(s)
        ret = c.fromdict(data[1])
      except Exception as e: 
        head = n + "### Core parser ###" + n 
        exc += head + spyder.exception()
        e.__context__ = None      
        raise
    except Exception as e: 
      try: #read using File, also works with web links
        try:
          ret = spyder.__types__["File"](filename, c).data()
        except Exception as e:
          head = n + "### Spyder.File parser ###" + n 
          exc += head + spyder.exception()
          e.__context__ = None
          raise                  
      except Exception as e:
        raise type(e)(exc)
    return ret
""" 
_tofilecode = """  def tofile(self, filename):
    \"\"\"Auto-generated by Spyder:
     module core
     file parse.py
     function generate_parse      
    Saves object to file\"\"\"
    f = spyder.__types__["File"](filename, type(self), "w", self)
    f.save()
    f.close()
"""

_totempfilecode = """  def totempfile(self):
    \"\"\"Auto-generated by Spyder:
     module core
     file parse.py
     function generate_parse      
    Saves object to temporary file\"\"\"
    import random, os,tempfile
    tempdir = tempfile.gettempdir()
    while 1:
      tempnam = tempdir + os.sep + str(random.randrange(1,1000000000)) +".tmp"
      if not spyder.loader.file_exists(tempnam): break
    f = spyder.__types__["File"](tempnam, type(self), "w", self)
    f.temporary(True)
    f.save()
    return f
"""      

def generate_parse(typename, parentnames, source, members, deleted_members, block):
  if block != None: raise Exception
  s = """  def __parse__(self, s):
    \"\"\"Auto-generated by Spyder:
     module core
     file parse.py
     function generate_parse
    Parsing constructor, which redirects to spyder.core.parse\"\"\"
    p = spyder.core.parse(s, \"Spyder\", True)
    if type(p) == type(self):
      self.__dict__ = p.__dict__
      return None
    if isinstance(p, Object):
      a = self.__spydercopy__(p)
    else:
      assert isinstance(p, tuple) and len(p) == 2, type(p) #if not a Spyder object, spyder.core.parse should return a tuple 
      args, kwargs = p[0], p[1]
      if len(kwargs) > 0 and len(args) == 0:
        a = self.__copydict__(p[1])
      else:
        a = self.__unpack__(*p[0], **p[1])
    return a
"""
  arraycode = """  def __parse__(self, s):
    \"\"\"Auto-generated by Spyder:
     module core
     file parse.py
     function generate_parse   
    Parsing constructor, which redirects to spyder.core.parse\"\"\"
    args = spyder.core.parse(s, \"Spyder\", True)
    self.__spydercopy__(args[0])
"""
  s += _fromfilecode
  s += _tofilecode
  s += _totempfilecode
  arraycode += _fromfilecode
  arraycode += _tofilecode
  arraycode += _totempfilecode
  if len(members) == 0: return None,arraycode    
  return s,arraycode

spyder.defineunimethod("__parse__", generate_parse)

