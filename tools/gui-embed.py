#!/usr/bin/env python2

"""
Command line tool to embed a Spyder model
"""

import sys, os
currdir = os.path.abspath(os.path.split(__file__)[0])
guidir = currdir + "/../gui"
sys.path.insert(0, guidir)
import spyder, Spyder
import attractmodel
from spyder.formtools import embed

inputfile = sys.argv[1]
inputdata = open(inputfile).read()
inputdir = os.path.split(inputfile)[0]
if len(inputdir): os.chdir(inputdir)
spydertypename, spyderdict = spyder.core.parse(inputdata)
spydertype = getattr(Spyder, spydertypename)
model = spydertype.fromdict(spyderdict)
embed(model)
print repr(model)