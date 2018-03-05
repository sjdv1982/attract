#!/usr/bin/env python

from serverconfig import *

def serve_upload():
  webform = cgi.FieldStorage()
  if "protocolfile" not in webform: raise Exception ### TODO, nice error message
  data = webform["protocolfile"].value
  interface = webform["interface"].value
  typ, content = spyder.core.parse(data)

  conversion = False
  if interface == "auto":
    interface = typ
  else:
    conversion = True
  if interface == "AttractModel":
    spydertype = Spyder.AttractModel
    formlib = form_model
    cgiscript = "attractserver.py"
  elif interface == "AttractEasyModel":
    spydertype = Spyder.AttractEasyModel
    formlib = form_standard
    cgiscript = "attractserver-easy.py"
  elif interface == "AttractPeptideModel":
    spydertype = Spyder.AttractPeptideModel
    formlib = form_peptide
    cgiscript = "attractserver-peptide.py"
  elif interface == "CryoInterface":
    spydertype = Spyder.CryoInterface
    formlib = form_cryo
    cgiscript = "attractserver-cryo.py"
  elif interface == "CryoInterface":
    spydertype = Spyder.CryoInterface
    formlib = form_cryo
    cgiscript = "attractserver-cryo.py"
  elif interface == "CryoEasyInterface":
    spydertype = Spyder.CryoEasyInterface
    formlib = form_cryo_easy
    cgiscript = "attractserver-cryo-easy.py"
  else:
    raise ValueError(interface) ### TODO, nice error message

  if conversion:
    spydertype2 = getattr(Spyder, typ)
    model = spydertype2.fromdict(content)
    model = model.convert(spydertype)
  else:
    model = spydertype.fromdict(content)
  import random
  from spyder.formtools import embed
  embed(model)
  mydir = "attract-upload" + str(random.randint(1,1000000))
  fname = "attract.web"
  os.chdir("/tmp")
  os.mkdir(mydir)
  os.chdir(mydir)
  model.tofile(fname)

  f = spydertype._form()
  f = formlib.webform(f, model)
  header = formlib.header.replace("<head>", "<head>\n        <base href=\"%s\" target=\"_blank\">" % webdir, 1)
  html = attracthtmlform.htmlform(
    obj=model,
    form=f,
    cgi=cgidir+cgiscript,
    header=header,
    footer=formlib.footer,
    header_indentation = 12,
    resourcefilename=mydir+"/"+fname,
  )
  return html

print "Content-type: text/html\n\n"

try:

  import os, sys, time, traceback, cgi
  import spyder, Spyder
  import spyder.formtools
  import attractmodel
  import attracthtmlform
  import form_model, form_standard, form_peptide, form_cryo, form_cryo_easy

  r = serve_upload()
  print r
except Exception, e:
  print "<pre>"
  print traceback.format_exc() + "\n"
  print "</pre>"
  sys.exit()
