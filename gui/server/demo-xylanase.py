#!/usr/bin/env python

from serverconfig import *

def serve_demo():
  filename = "demo/demo-xylanase.web"
  assert os.path.exists(filename)

  spydertype = Spyder.AttractEasyModel    
  formlib = form_standard
  cgiscript = "attractserver-easy.py"
  model = spydertype.fromfile(filename)
  
  import random
  from spyder.formtools import embed
  embed(model)     
  filename2 = os.path.split(filename)[0] + "-" + str(random.randint(1,1000000))  
  os.chdir("/tmp")
  model.tofile(filename2)
  
  f = spydertype._form()
  f = formlib.webform_easy(f, model)
  header = formlib.header.replace("<head>", "<head>\n        <base href=\"%s\" target=\"_blank\">" % webdir, 1)
  html = attracthtmlform.htmlform(
    obj=model, 
    form=f, 
    cgi=cgidir+cgiscript, 
    header=header, 
    footer=formlib.footer, 
    header_indentation = 12,
    resourcefilename=filename2,
  )
  return html

print "Content-type: text/html\n\n"

try:
  
  import os, sys, time, traceback, cgi
  import spyder, Spyder
  import spyder.formtools
  import attractmodel
  import attracthtmlform
  import form_model, form_standard

  r = serve_demo()
  print r
except Exception, e:  
  print "<pre>"
  print traceback.format_exc() + "\n"
  print "</pre>"  
  sys.exit()


