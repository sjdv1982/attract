#!/usr/bin/env python
  
print "Content-type: text/html\n\n<pre>"

try:
  
  import os, sys, time, traceback
  sys.path.append(os.path.split(os.path.abspath(__file__))[0] + "/..") 
  import spyder, Spyder
  import spyder.cgimodule
  import attractmodel
  import formeasy
  import attractsave
  from serverlib import serve_attract

  r = serve_attract(Spyder.AttractEasyModel, formeasy, attractsave.deploy_easy)
  print r
except Exception, e:
  print traceback.format_exc() + "\n"
  print "</pre>"  
  sys.exit()

print "</pre>"  
