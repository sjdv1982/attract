#!/usr/bin/env python3

import os
import sys
currdir = os.path.abspath(os.path.split(__file__)[0])
guidir = currdir + "/../gui"
# Ensure local SpyderSilk is importable (provides `spyder` and registers `Spyder`)
spydersilkdir = os.path.abspath(os.path.join(currdir, "..", "spydersilk-0.03"))
sys.path.insert(0, spydersilkdir)
sys.path.insert(0, guidir)

import spyder  # noqa: E402
import Spyder  # noqa: E402
from spyder.formtools import make_relpath, check_embedded
import attractmodel
import attractsave

inputfile = sys.argv[1]
inputfile_abs = os.path.abspath(inputfile)
inpdir = os.path.abspath(os.path.split(inputfile_abs)[0])
os.chdir(inpdir)
with open(inputfile_abs, "r", encoding="utf-8") as f:
  spydertypename, data = spyder.core.parse(f.read())
spydertype = getattr(Spyder, spydertypename)
model = spydertype.fromdict(data)
if spydertypename == "AttractPeptideModel": 
  model = model.convert(Spyder.AttractEasyModel)  
  r = model.partners[1].pdbfile
  r.link("./peptide.pdb")
  r.save()
  
embedded = check_embedded(model)
if embedded is not None:
  # If the model embeds resources (common for demos), deploy them next to the
  # input file so we can generate a runnable script.
  if spydertypename == "AttractModel":
    deploy = attractsave.deploy
  elif spydertypename == "AttractEasyModel":
    deploy = attractsave.deploy_easy
  elif spydertypename == "AttractPeptideModel":
    deploy = attractsave.deploy_peptide
  else:
    deploy = None

  if deploy is None:
    print("Cannot generate shell script: %s is an embedded resource, not a file name" % embedded)
    sys.exit(1)

  deploy(model, inpdir)
  embedded = check_embedded(model)
  if embedded is not None:
    print("Cannot generate shell script: %s is an embedded resource, not a file name" % embedded)
    sys.exit(1)
make_relpath(inpdir, model)
sh = os.path.splitext(inputfile_abs)[0] + ".sh"
script = model.generate()
fsh = open(sh, "w")
fsh.write(script+"\n")
fsh.close()
os.system("chmod +x %s" % sh)
print("Shell script generated: %s" % sh)
