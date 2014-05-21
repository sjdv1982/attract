webdir = "http://www.attract.ph.tum.de/services/ATTRACT/"
localdir = "/home/server/services/ATTRACT/html/"
resultdir = "results/"

import os, sys, cgi, datetime

def serve_attract(spydertype, formlib, deploy):
  import spyder
  import spyder.htmlform
  from spyder.formtools import embed
  import pprint
  import traceback
  webform = spyder.htmlform.dict_from_fieldstorage(cgi.FieldStorage())
  f = formlib.webserverform(webform, spydertype=spydertype)
  
  resourcemodel = None  
  if "_tmpresource" in webform:
    tmpf = "/tmp/" + webform["_tmpresource"]
    resourcemodel = spydertype.fromfile(tmpf)
  webdict = spyder.htmlform.cgi(webform,f,resourcemodel)  
  try:
    newmodel = spydertype.fromdict(webdict)
  except:    
    raise ValueError(traceback.format_exc() + "\n" + pprint.pformat(webform) + "\n" + pprint.pformat(webdict))
  runname = newmodel.convert(Spyder.AttractModel).runname
  if runname is None: runname = (datetime.datetime.now()).isoformat()
  import string
  runname = ''.join([str(char) for char in runname if char in string.printable])
  import re
  runname = re.sub(r'[^a-zA-Z0-9\._-]', '', runname)
  newmodel.runname = runname
  embed(newmodel)
  import random
  mydir = "run" + str(random.randint(1,1000000))  
  fname = "attract-%s.web" % runname
  os.chdir(localdir + resultdir)
  os.mkdir(mydir)
  os.chdir(mydir)
  newmodel.tofile(fname)
  os.mkdir("attract-%s" %runname)
  os.chdir("attract-%s" %runname)
  deploy(newmodel, ".")
  newmodel.tofile(fname)
  script = newmodel.generate()
  f = open("attract-%s.sh" %runname, "w")
  f.write(script)
  f.close()
  os.system("chmod +x attract-%s.sh" %runname)
  os.chdir("..")
  ret = []
  ret.append("You can download your self-contained parameter file <a href='%s'>here</a>" % (webdir+resultdir+mydir+"/"+fname))
  aname = "attract-%s.tgz" % runname
  os.system("tar czf %s attract-%s" % aname %runname)
  ret.append("You can download an archive with executable shell script, all docking files, and linked parameter file <a href='%s'>here</a>" % (webdir+resultdir+mydir+"/"+aname))
  return "\n".join(ret)
