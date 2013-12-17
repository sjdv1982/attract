#!/usr/bin/env python

webdir = "http://www.attract.ph.tum.de/services/ATTRACT/"
localdir = "/home/server/services/ATTRACT/html/"
resultdir = "results/"

def serve_attract():
  webform = spyder.htmlform.dict_from_fieldstorage(cgi.FieldStorage())
  f = form_model.webserverform(webform, spydertype=Spyder.AttractModel)
  
  resourcemodel = None  
  if "_tmpresource" in webform:
    tmpf = webform["_tmpresource"].value    
    resourcemodel = Spyder.AttractModel.fromfile(tmpf)    
    nr_iterations = max( len(resourcemodel.iterations) , nr_iterations )    
  webdict = spyder.htmlform.cgi(webform,f,resourcemodel)  
  newmodel = Spyder.AttractModel.fromdict(webdict)
  resources.embed(newmodel)
  import random
  mydir = "run" + str(random.randint(1,1000000))
  fname = "attract.web"
  os.chdir(localdir + resultdir)
  os.mkdir(mydir)
  os.chdir(mydir)
  newmodel.tofile(fname)
  os.mkdir("attract")
  os.chdir("attract")
  resources.deploy(newmodel, ".")
  newmodel.tofile(fname)
  script = newmodel.generate()
  f = open("attract.sh", "w")
  f.write(script)
  f.close()
  os.system("chmod +x attract.sh")
  os.chdir("..")
  ret = []
  ret.append("You can download your self-contained parameter file <a href='%s'>here</a>" % (webdir+resultdir+mydir+"/"+fname))
  aname = "attract.tgz"
  os.system("tar czf %s attract" % aname)
  ret.append("You can download an archive with executable shell script, all docking files, and linked parameter file <a href='%s'>here</a>" % (webdir+resultdir+mydir+"/"+aname))
  return "\n".join(ret)
  
print "Content-type: text/html\n\n<pre>"

try:
  
  import os, sys, time, traceback, cgi
  import spyder, Spyder
  import spyder.htmlform
  import attractmodel
  import form_model
  import resources

  r = serve_attract()
  print r
except Exception, e:
  print traceback.format_exc() + "\n"
  print "</pre>"  
  sys.exit()

print "</pre>"  
