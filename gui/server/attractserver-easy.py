#!/usr/bin/env python
from __future__ import print_function

print("Content-type: text/html\n\n<pre>")

try:  
  import os, sys, time, traceback
  sys.path.append(os.path.split(os.path.abspath(__file__))[0] + "/..") 
  import spyder, Spyder
  import spyder.formtools
  import attractmodel
  import formeasy
  import attractsave
  from serverlib import serve_attract, AttractServerError
  os.system("chmod a+r+w /tmp/*spy.*")
except Exception, e:
  s = "There was an error in the initialization of the server\n"
  s += """Please save this page from your browser as HTML and email it to sjoerd@tum.de, and it will be fixed as soon as possible.
The ATTRACT web interface is in active development, thank you for your patience.  
"""
  s += "\n<B>Full error information</B>:\n\n"
  s += traceback.format_exc()
  print(s)
  print("</pre>")
  sys.exit()
  
try:  
  r = serve_attract(Spyder.AttractEasyModel, formeasy, attractsave.deploy_easy)
  print(r)
except AttractServerError as e:
  print("There was an inconsistency in your data")
  print("\n<B>Error message</B>")
  print(e.status)
  if e.delta is not None:
    print("\n<B>You can download the following file</B>")
    print("\n")
    print("\n".join(e.delta))
  print("</pre>") 
  sys.exit()  
except Exception as e:
  print("There was an unknown error in the server")
  print("""Please save this page from your browser as HTML and email it to sjoerd@tum.de, and it will be fixed as soon as possible.
The ATTRACT web interface is in active development, thank you for your patience.  
""")  
  print(traceback.format_exc() + "\n")
  print("</pre>") 
  sys.exit()

print("</pre>")
