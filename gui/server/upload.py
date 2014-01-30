#!/usr/bin/env python

webdir = "http://www.attract.ph.tum.de/services/ATTRACT/"
cgidir = "http://www.attract.ph.tum.de/cgi/services/ATTRACT/"
localdir = "/home/server/services/ATTRACT/html/"
resultdir = "results/"
cgiscript = "attractserver.py"

def serve_upload():
  webform = cgi.FieldStorage()
  if "protocolfile" not in webform: raise Exception ### TODO, nice error message
  data = webform["protocolfile"].value
  typ, content = spyder.core.parse(data)
  if typ != "AttractModel": raise ValueError(typ) ### TODO, nice error message
  model = attractmodel.AttractModel.fromdict(content)
  import random
  from spyder.formtools import embed
  embed(model)   
  mydir = "run" + str(random.randint(1,1000000))  
  fname = "attract.web"
  os.chdir(localdir + resultdir)
  os.mkdir(mydir)
  os.chdir(mydir)
  model.tofile(fname)
  
  for p in model.partners: 
    p.rmsd, p.use_modes = None, None #"virtual" form attributes defined in form.webform
  f = attractmodel.AttractModel._form()
  f = form.webform(f, model)
  header = form.header.replace("<head>", "<head>\n        <base href=\"%s\" target=\"_blank\">" % webdir, 1)
  html = attracthtmlform.htmlform(
    obj=model, 
    form=f, 
    cgi=cgidir+cgiscript, 
    header=header, 
    footer=form.footer, 
    header_indentation = 12,
    hidden = {"_tmpresource":localdir+resultdir+mydir+"/"+fname},
  )
  return html
  
  #TODO: generate HTML instead with "_tmpresource" attribute
  #return "You can download your parameter file <a href='%s'>here</a>" % (webdir+resultdir+mydir+"/"+fname)

print "Content-type: text/html\n\n"

try:
  
  import os, sys, time, traceback, cgi
  import spyder, Spyder
  import spyder.htmlform
  import attractmodel
  import attracthtmlform
  import form

  r = serve_upload()
  print r
except Exception, e:  
  print "<pre>"
  print traceback.format_exc() + "\n"
  print "</pre>"  
  sys.exit()


