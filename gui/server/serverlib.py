from serverconfig import *
import os, sys, cgi, datetime, string, re, random
import spyder
import spyder.formtools
from spyder.formtools import embed
import pprint
import traceback

class AttractServerError(Exception):
  def __init__(self, status, delta):
    Exception.__init__(self)
    self.status = status
    self.delta = delta

def format_delta(delta):
  if delta is None: return None
  return ["DELTA:\n" + str(delta)]
  
def serve_attract(spydertype, formlib, deploy):
  webdict = spyder.formtools.cgi.dict_from_fieldstorage(cgi.FieldStorage())
  f = formlib.webserverform(webdict, spydertype=spydertype)
  
  resourceobj = None  
  resourcefilevar = getattr(f, "resourcefilevar", None)
  if resourcefilevar is not None and resourcefilevar in webdict:
    tmpf = "/tmp/" + webdict[resourcefilevar]
    resourceobj = spydertype.fromfile(tmpf)
  newmodel, status, delta = spyder.formtools.cgi.cgi(webdict, f, resourceobj, spydertype=spydertype)      
 
  cwd = os.getcwd()
  os.chdir(localresultdir)  
  mydir = "run" + str(random.randint(1,1000000))    
  os.mkdir(mydir)
  
  if delta is None:
    raise AttractServerError(status="You didn't submit any data", delta=None)
      
  if newmodel is None or status.splitlines()[0].find("OK") == -1:    
    raise AttractServerError(status=status, delta=format_delta(delta))
  
  os.chdir(cwd)
  runname = getattr(newmodel, "runname", None)
  if runname is None or runname == "attract": runname = "attract-" + datetime.datetime.now().date().isoformat()
  runname = ''.join([str(char) for char in runname if char in string.printable])
  runname = re.sub(r'[^a-zA-Z0-9\.]', '_', runname)
  newmodel.runname = runname
  embed(newmodel)
  fname_embedded = "%s-embedded.web" % runname
  fname = "%s.web" % runname  
  os.chdir(localresultdir)
  os.chdir(mydir)
  deltafile = format_delta(runname)
  newmodel.tofile(fname_embedded)
  os.mkdir(runname)
  os.chdir(runname)
  deploy(newmodel, ".")
  newmodel.tofile(fname)
  script = newmodel.generate()
  f = open("%s.sh" % runname, "w")
  f.write(script)
  f.close()
  os.system("chmod +x %s.sh" %runname)
  os.chdir("..")
  os.system("chmod a+rw %s" % mydir)
  ret = []
  ret += format_delta(delta)
  ret += ["",status,""]
  ret += ["",str(newmodel),""]
  ret.append("You can download your self-contained parameter file <a href='%s'>here</a>" % (webresultdir+mydir+"/"+fname_embedded))
  aname = "%s.tgz" % runname
  os.system("tar czf %s %s" % (aname , runname))
  ret.append("You can download an archive with executable shell script, all docking files, and linked parameter file <a href='%s'>here</a>" % (webresultdir+mydir+"/"+aname))
  return "\n".join(ret)
