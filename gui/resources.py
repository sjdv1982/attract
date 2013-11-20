import Spyder

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

def embed(m):
  try:
    items = [a for a in m]
    items = enumerate(items)
    ar = True
  except TypeError:
    items = list(m.__dict__.items())
    ar = False
  for k,v in items:
    if isinstance(v, Spyder.Resource): 
      v.embed()
    elif isinstance(v, Spyder.Object):
     embed(v)      
