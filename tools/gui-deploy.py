#!/usr/bin/env python2

"""
Command line tool to deploy an embedded ATTRACT Spyder model into a directory
usage: gui-deploy.py <inputfile> <directory>
"""

import sys, os
currdir = os.path.abspath(os.path.split(__file__)[0])
guidir = currdir + "/../gui"
sys.path.insert(0, guidir)
import spyder, Spyder
from spyder.formtools import make_relpath
import attractmodel
import attractsave

inputfile = sys.argv[1]
directory = sys.argv[2]
assert os.path.isdir(directory)
directory = os.path.abspath(directory)

inputdata = open(inputfile).read()
inputdir = os.path.split(inputfile)[0]
if len(inputdir): os.chdir(inputdir)
spydertypename, spyderdict = spyder.core.parse(inputdata)
spydertype = getattr(Spyder, spydertypename)
if spydertypename == "AttractModel":
  deploy = attractsave.deploy
elif spydertypename == "AttractEasyModel":
  deploy = attractsave.deploy_easy
elif spydertypename == "AttractPeptideModel":
  deploy = attractsave.deploy_peptide
else:
  raise Exception("Don't know how to deploy Spyder model '%s'" % spydertypename)
model = spydertype.fromdict(spyderdict)
deploy(model, directory)
make_relpath(inputdir, model)
print repr(model)
