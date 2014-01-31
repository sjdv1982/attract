import os
import Spyder

from spyder.formtools import make_relpath, check_embedded
  
def save(m, outp, *args):
  v = m._get()   
  outpdir = os.path.split(outp)[0]  
  if v is not None:
    make_relpath(outpdir, v)
    v.tofile(outp)
    print("Form object saved to %s" % outp)
  else:
    print("Cannot save form: does not contain a valid object")

def generate(m, outp, *args):
  v = m._get() 
  outpdir = os.path.split(outp)[0]
  if v is not None:
    embedded = check_embedded(v)
    if embedded is not None:
      print("Cannot generate shell script: %s is an embedded resource, not a file name" % embedded)
      return
    make_relpath(outpdir, v)
    sh = os.path.splitext(outp)[0] + ".sh"
    script = v.generate()
    fsh = open(sh, "w")
    fsh.write(script+"\n")
    fsh.close()
    os.system("chmod +x %s" % sh)
    print("Shell script generated: %s" % sh)
  else:
    print("Cannot generate shell script: form does not contain a valid object")
    
def _deploy(resource, fname):
  if resource is not None: 
    resource.link(fname)
    resource.save()
  
def deploy(model, dir):
  d = dir + "/"
  if dir in (None, "", ".", "./"): d = ""
  elif dir.endswith("/"): d = dir
  for n,p in enumerate(model.partners):
    _deploy(p.pdbfile,d+"partner-%d.pdb" % (n+1))
    _deploy(p.ensemble_list,d+"ensemble-%d.list" % (n+1))
    _deploy(p.modes_file,d+"partner-%d.modes" % (n+1))
    _deploy(p.aa_modes_file,d+"partner-aa-%d.modes" % (n+1))
    _deploy(p.rmsd_pdb,d+"partner-rmsd-%d.pdb" % (n+1))
    _deploy(p.collect_pdb,d+"partner-collect-%d.pdb" % (n+1))
    _deploy(p.collect_ensemble_list,d+"partner-collect-ensemble-%d.list" % (n+1))
  _deploy(model.cryoem_data,d+"cryo.map")
  _deploy(model.cryoem_scoring_data,d+"cryo-scoring.map")
  _deploy(model.start_structures_file,d+"startstruc.dat")
  _deploy(model.rotations_file,d+"rotations.dat")
  _deploy(model.translations_file,d+"translations.dat")

def deploy_easy(model, dir):
  d = dir + "/"
  if dir in (None, "", ".", "./"): d = ""
  elif dir.endswith("/"): d = dir
  for n,p in enumerate(model.partners):
    _deploy(p.pdbfile,d+"partner-%d.pdb" % (n+1))
    _deploy(p.rmsd_pdb,d+"partner-rmsd-%d.pdb" % (n+1))
  