#!/usr/bin/env python

webdir = "http://localhost/services/attractcgi/"
cgidir = "http://localhost/cgi/services/attractcgi/"
localdir = "/home/sjoerd/services/attractcgi/html/"
resultdir = "results/"
cgiscript = "attractserver.py"

def serve_upload():
  webform = cgi.FieldStorage()
  if "protocolfile" not in webform: raise Exception ### TODO, nice error message
  data = webform["protocolfile"].value
  typ, content = spyder.core.parse(data)
  if typ != "AttractModel": raise ValueError(typ) ### TODO, nice error message
  model = attractmodel.AttractModel.fromdict(content)
  resources.embed(model)
  import random
  mydir = "run" + str(random.randint(1,1000000))  
  fname = "attract.web"
  os.chdir(localdir + resultdir)
  os.mkdir(mydir)
  os.chdir(mydir)
  model.tofile(fname)
  
  for p in model.partners: p.rmsd = None #"virtual" form attribute defined in form.webform
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
  import resources

  r = serve_upload()
  print r
except Exception, e:  
  print "<pre>"
  print traceback.format_exc() + "\n"
  print "</pre>"  
  sys.exit()


