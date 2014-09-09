from serverconfig import *
import os, sys, cgi, datetime, string, re, random
import spyder
import spyder.formtools
from spyder.formtools import embed
import pprint, json
import traceback

response_success = """<B>Your docking parameters were received successfully</B>
A docking protocol script was generated from your parameters

<!--
You will now receive three files for download:
//-->
"""

response_processing1_error = """<B>Your docking parameters were received successfully, but they could not be processed</B>
This is probably a bug in the ATTRACT server.
Please send the error message below and your delta file to sjoerd@tum.de, and it will be fixed as soon as possible

The ATTRACT web interface is in active development, thank you for your patience.

<B>Error message</B>
%s

"""

response_processing2_error = """<B>Your docking parameters were received successfully, but they could not be processed</B>
This is probably a bug in the ATTRACT server.
Please send the error message below, your parameter file and your delta file to sjoerd@tum.de, and it will be fixed as soon as possible

The ATTRACT web interface is in active development, thank you for your patience.

<B>Error message</B>
%s
<!--
You will now receive two files for download:
//-->
"""

response_generator_error = """<B>Your docking parameters were received successfully, but no protocol script could be generated</B>
This is probably a bug in the ATTRACT script generator. 
Please send the error message below, your parameter file and your delta file to sjoerd@tum.de, and it will be fixed as soon as possible

The ATTRACT web interface is in active development, thank you for your patience.  

<B>Error message</B>
%s

<!--
You will now receive two files for download:
//-->

"""

response_directory =  """<B>Docking directory</B>
The docking directory is provided as a tar-zipped archive (.tgz format)
The directory contains a shell script (<B>%s</B>) that you can execute directly on your computer.
This performs the docking and analysis using ATTRACT

<b><i>Download the docking directory: <a href='%s'>%s</a></i></b>

<!--
The directory also contains a deployed parameter file (%s) that can be edited locally. 
Its resources (PDB files etc.) are deployed: they refer to the file names in the docking directory.

<u>Deployed parameter files in the local ATTRACT GUI</u>
You can edit the deployed parameter file directly with the local ATTRACT GUI, and generate a shell script

<u>Advanced usage: deployed parameter files in the web GUI</u>
You cannot upload deployed parameter files to the web interface. You have to embed your parameter file using gui-embed, or use the embedded parameter file instead.
//-->

"""

response_embedded = """<!-- <B>Embedded parameter file</B>
The embedded parameter file contains all docking parameters needed by ATTRACT, in a single file
This file describes your docking protocol in a reproducible manner, e.g. for automatic recalculation

<u>Embedded parameter files in the web GUI</u>
You can upload the embedded parameter file to the web interface and modify its parameters

Download the embedded parameter file: <a href='%s'>%s</a>

<u>Advanced usage: embedded parameter files in the local ATTRACT GUI</u>
You can edit the embedded parameter file directly with the local ATTRACT GUI
However, this file contains embedded resources (PDB files etc.). If you want to generate a new docking script, you have to deploy it first into a directory using gui-deploy, or use the deployed parameter file instead
//-->

"""

response_deltafile = """<!-- <B>Delta file</B>
The delta file contains the web form parameters that were filled in or changed. This file is in text format (JSON).
The delta file is the most useful reference file for describing your docking protocol in words, e.g. a Materials & Methods section
Providing your delta file is essential for help, support, feedback and bug reports.

Download the delta file: <a href='%s'>%s</a>

<u>Advanced usage: delta files in the local ATTRACT GUI</u>
You can load a delta file with the local ATTRACT GUI with the option: --delta &ltdelta file&gt
Delta files contain embedded resources (PDB files etc.): if you want to generate a new docking script, you have to save it as a parameter file and deploy it first into a directory using gui-deploy
//-->
"""

class AttractServerError(Exception):
  def __init__(self, status, delta):
    Exception.__init__(self)
    self.status = status
    self.delta = delta

def format_delta(delta, mydir, runname="attract"):
  if delta is None: return ""
  fname = runname + "-delta.json"
  try:
    outputf = open(mydir+"/"+fname, "w") 
    json.dump(delta,outputf,sort_keys=True,indent=4)
    outputf.close()  
    return response_deltafile % (webresultdir+mydir+"/"+fname, fname)
  except:
    return "<B>Delta</B>\n" + pprint.pformat(delta)
  
  
def serve_attract(spydertype, formlib, deploy):
  # Obtain the webdict (what the user submitted) and f, the spyderform used to generate the HTML
  webdict = spyder.formtools.cgi.dict_from_fieldstorage(cgi.FieldStorage())
  f = formlib.webserverform(webdict, spydertype=spydertype)
  
  # Fill in links to embedded resources
  resourceobj = None  
  resourcefilevar = getattr(f, "resourcefilevar", None)
  if resourcefilevar is not None and resourcefilevar in webdict:
    tmpf = "/tmp/" + webdict[resourcefilevar]
    resourceobj = spydertype.fromfile(tmpf)
  newmodel, status, delta = spyder.formtools.cgi.cgi(webdict, f, resourceobj, spydertype=spydertype)      
 
  # Detect empty form  
  if not len(webdict) or delta is None:
    raise AttractServerError(status="You didn't submit any data", delta=None)      
  
  # Create a result directory
  cwd = os.getcwd()
  os.chdir(localresultdir)  
  mydir = "run" + str(random.randint(1,1000000))    
  os.mkdir(mydir)
  
  # Detect errors in the form submission
  # We will still give the user the delta file 
  deltamessage = format_delta(delta, mydir)
  if status.splitlines()[0].find("OK") == -1:    
    raise AttractServerError(status=status, delta=deltamessage)
  
  # New model was received OK, should be no bugs from here on...  
  try:
    # Obtain a runname
    runname = getattr(newmodel, "runname", None)
    if runname is None or runname == "attract": runname = "attract-" + datetime.datetime.now().date().isoformat()
    runname = ''.join([str(char) for char in runname if char in string.printable])
    runname = re.sub(r'[^a-zA-Z0-9\.]', '_', runname)
    
    # Embed the model, and save the delta
    os.chdir(cwd)
    newmodel.runname = runname
    embed(newmodel)  
    fname_embedded = "%s-embedded.web" % runname
    fname = "%s.web" % runname  
    os.chdir(localresultdir)    
    deltamessage = format_delta(delta, mydir, runname)
    os.chdir(mydir)
    newmodel.tofile(fname_embedded)
  except: 
    # Processing error 1
    error = traceback.format_exc()
    response = response_processing1_error % error
    response += deltamessage
    return response
    
  try:  
    #Create an archive directory equal to runname
    os.mkdir(runname)
    os.chdir(runname)
    deploy(newmodel, ".")
    newmodel.tofile(fname)
  except: 
    # Processing error 2
    error = traceback.format_exc()
    response = response_processing2_error % error
    response += response_embedded % (webresultdir+mydir+"/"+fname_embedded, fname_embedded)
    response += deltamessage
    return response
  
  try:
    #Run the generator, write the script and create the archive
    script = newmodel.generate()
    shname = "%s.sh" % runname
    f = open(shname, "w")
    f.write(script)
    f.close()
    os.system("chmod +x %s" % shname)
    os.chdir("..")
    os.system("chmod a+rw %s" % mydir)
    aname = "%s.tgz" % runname
    os.system("tar czf %s %s" % (aname,runname))    
    assert os.path.exists(aname)
  except:
    # Generator error
    error = traceback.format_exc()
    response = response_generator_error % error
    response += response_embedded % (webresultdir+mydir+"/"+fname_embedded, fname_embedded)
    response += deltamessage
    return response

  response = response_success
  response += response_directory % (shname, webresultdir+mydir+"/"+aname, aname, fname)
  response += response_embedded % (webresultdir+mydir+"/"+fname_embedded, fname_embedded)
  response += deltamessage
  return response
    