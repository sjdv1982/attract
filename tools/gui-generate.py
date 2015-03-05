import os, sys
currdir = os.path.abspath(os.path.split(__file__)[0])
guidir = currdir + "/../gui"
sys.path.insert(0, guidir)
import spyder, Spyder
from spyder.formtools import make_relpath, check_embedded
import attractmodel
from attractsave import deploy_peptide

inputfile = sys.argv[1]
inpdir = os.path.abspath(os.path.split(inputfile)[0])
os.chdir(inpdir)
spydertypename, data = spyder.core.parse(open(inputfile).read())
spydertype = getattr(Spyder, spydertypename)
model = spydertype.fromdict(data)
if spydertypename == "AttractPeptideModel": 
  model = model.convert(Spyder.AttractEasyModel)  
  r = model.partners[1].pdbfile
  r.link("./peptide.pdb")
  r.save()
  
embedded = check_embedded(model)
if embedded is not None:
  print("Cannot generate shell script: %s is an embedded resource, not a file name" % embedded)
  sys.exit()
make_relpath(inpdir, model)
sh = os.path.splitext(inputfile)[0] + ".sh"
script = model.generate()
fsh = open(sh, "w")
fsh.write(script+"\n")
fsh.close()
os.system("chmod +x %s" % sh)
print("Shell script generated: %s" % sh)



