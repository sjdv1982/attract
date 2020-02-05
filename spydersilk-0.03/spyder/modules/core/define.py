# Copyright 2007-2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 
  
def statement_define(l, title, block):
  ret = '"""Define %s"""\n' % title.replace('"""', "\\\"\\\"\\\"")
  if block != None: ret = '"""Define %s {' % title + block.replace('"""', "\\\"\\\"\\\"") + '}"""\n'
  supported = []
  start = -1
  end = -1
  for n in range(len(title)):
    if start == -1 and title[n] == "(": start = n
    if end == -1 and title[len(title)-n-1] == ")": end = len(title)-n-1
  before = title[:start].strip()
  if len(before.split()) > 1:
    raise Exception("compile error: malformed Define statement: %s" % title)
  between = title[start+1:end].split(",")
  for n in range(len(between)):
    between[n] = between[n].split()
    for nn in range(len(between[n])): between[n][nn] = between[n][nn].strip()
  after = title[end+1:].strip()
  if (block == None) != (len(after) > 0):
    raise Exception("compile error: malformed Define statement: %s" % title)
  if (len(between[0]) == 1) != (block ==None) or (len(between[0]) > 2):
    raise Exception("compile error: malformed Define statement: %s" % title)
  if block != None:
    spaceblock = ""
    spaces = -1
    for l in block.split('\n'):
      if len(l.strip()) == 0: continue
      if spaces == -1:
        spaces = len(l) - len(l.lstrip())
      spaceblock += "  %s\n" % l[spaces:]    
    if spaces == -1: spaceblock += "    pass\n"
  if len(between) == 1:
    if len(between[0]) == 1:      
      if after in ("CAST", "SPLIT"): after = "\"" + after + "\""
      else:
        after = after.split()
        if len(after) != 1 and (not after[1].startswith('"""') or not after[-1].rstrip().endswith('"""')):
          raise Exception("compile error: malformed Define statement: %s" % title)
        after = after[0]
      ret += "spyder.core.defineconverter(\"%s\",\"%s\",%s)\n" % (between[0][0], before,after)
    else:
      ret += "def spyderconverterfunction_%d(%s):\n%s\n" % (spyder.core.define_functioncounter, between[0][1],spaceblock)
      ret += "spyder.core.defineconverter(\"%s\",\"%s\",spyderconverterfunction_%d)\n\ndel spyderconverterfunction_%d\n" % (between[0][0], before,spyder.core.define_functioncounter,spyder.core.define_functioncounter)
    supported.append(before)
    supported.append(between[0][0])
  else:
    raise Exception("compile error: malformed Define statement: %s" % title)    
  if block != None: spyder.core.define_functioncounter += 1
  spyder.safe_eval(ret)
  return ret,supported

def statement_method(l, title, block):
  ret = '"""Method %s\n"""\n' % title.replace('"""', "\\\"\\\"\\\"")
  if block != None: ret = '"""Method %s {' % title + block.replace('"""', "\\\"\\\"\\\"") + '}"""\n'
  supported = []
  start = -1
  end = -1
  for n in range(len(title)):
    if start == -1 and title[n] == "(": start = n
    if end == -1 and title[len(title)-n-1] == ")": end = len(title)-n-1
  before = title[:start].strip()
  if len(before.split()) > 1:
    raise Exception("compile error: malformed Method statement: %s" % title)
  between = title[start+1:end].split(",")
  for n in range(len(between)):
    between[n] = between[n].split()
    for nn in range(len(between[n])): between[n][nn] = between[n][nn].strip()
  after = title[end+1:].strip()
  if (block == None) != (len(after) > 0):
    raise Exception("compile error: malformed Method statement: %s" % title)
  if (len(between[0]) == 1) != (block == None) or (len(between[0]) > 2):
    raise Exception("compile error: malformed Method statement: %s" % title)
  if block != None:
    spaceblock = ""
    spaces = -1
    for l in block.split('\n'):
      if len(l.strip()) == 0: continue
      if spaces == -1:
        spaces = len(l) - len(l.lstrip())
      spaceblock += "  %s\n" % l[spaces:]    
    if spaces == -1: spaceblock += "    pass\n"
  if len(between) == 1:
    if len(between[0]) == 1:
      after = after.split()
      if len(after) != 1 and (not after[1].startswith('"""') or not after[-1].rstrip().endswith('"""')):
        raise Exception("compile error: malformed Method statement: %s" % title)
      after = after[0]
      ret += "spyder.core.definemethod(\"%s\",\"%s\",%s)\n" % (before, between[0][0], after)
    else:
      ret += "def spyderdefinemethod_%d(%s):\n%s\n" % (spyder.core.define_functioncounter, between[0][1],spaceblock)
      ret += "spyder.core.definemethod(\"%s\",\"%s\",spyderdefinemethod_%d)\ndel spyderdefinemethod_%d\n" % (before, between[0][0], spyder.core.define_functioncounter, spyder.core.define_functioncounter)
    supported.append(between[0][0])
  else:
    raise Exception("compile error: malformed Method statement: %s" % title)    
  if block != None: spyder.core.define_functioncounter += 1
  spyder.safe_eval(ret)
  return ret,supported

spyder.definestatement("Define", statement_define)
spyder.definestatement("Method", statement_method)
