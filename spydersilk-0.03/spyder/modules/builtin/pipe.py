# Copyright 2007-2011,2013, Sjoerd de Vries
# This file is part of the Spyder module: "builtin" 
# For licensing information, see LICENSE.txt 
from subprocess import Popen as _Popen, PIPE as _PIPE
def pipe(cmd, stdin=None):
  import  os, select
  bufsize = 3000
  p = _Popen(cmd,shell=True,bufsize=bufsize,
   stdin=_PIPE,stdout=_PIPE,stderr=_PIPE,close_fds=True)
  i,o,e = p.stdin,p.stdout,p.stderr
  stdout = ""
  stderr = ""
    
  ni,no,ne = i.fileno(), o.fileno(), e.fileno()  
  
  if stdin != None:
    size = 0
    for n in range(0,len(stdin),bufsize):
      try:
        os.write(ni, stdin[n:n+bufsize])
      except OSError:
        break
      ii = []
      while len(ii) == 0:        
        oo,ii,ee = select.select([o,e],[i],[])        
        change = False
        if o in oo: 
          stdout += os.read(no,bufsize)
          change = True
        if e in oo: 
          stderr += os.read(ne, bufsize)          
          change = True
        if change == False: break
  i.close()  
  oo,ii,ee = select.select([o,e],[],[])                        
  
  change = True  
  while (len(oo) > 0 or len(ee) > 0) and change:                      
    change = False
    stdoutlen = len(stdout)
    if o in oo: 
      stdout += os.read(no, bufsize)        
      if len(stdout) > stdoutlen: change = True
    stderrlen = len(stderr)
    if e in oo: 
      stderr += os.read(ne, bufsize)
      if len(stderr) > stderrlen: change = True
    oo,ii,ee = select.select([o,e],[],[])                

  return stdout,stderr
