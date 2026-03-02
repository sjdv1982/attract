#!/usr/bin/env python3

"""
Command line tool to deploy an embedded ATTRACT Spyder model into a directory
usage: gui-deploy.py <inputfile> <directory>
"""

import os
import sys
currdir = os.path.abspath(os.path.split(__file__)[0])
guidir = currdir + "/../gui"
# Ensure local SpyderSilk is importable (provides the `spyder` package and
# registers the dynamic `Spyder` module alias as a side effect).
spydersilkdir = os.path.abspath(os.path.join(currdir, "..", "spydersilk-0.03"))
sys.path.insert(0, spydersilkdir)
sys.path.insert(0, guidir)

import spyder  # noqa: E402
import Spyder  # noqa: E402
from spyder.formtools import make_relpath
import attractmodel
import attractsave

inputfile = sys.argv[1]
directory = sys.argv[2]
assert os.path.isdir(directory)
directory = os.path.abspath(directory)

with open(inputfile, "r", encoding="utf-8") as f:
  inputdata = f.read()
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
print(repr(model))
