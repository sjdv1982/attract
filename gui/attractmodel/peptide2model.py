peptide2model_version = "Converted from AttractPeptideModel by peptide2model 1.1"
import os
import sys

def peptide2model(pmodel):
  """
  Converts AttractPeptideModel to AttractEasyModel
  """
  partners = []
  #Convert receptor protein
  p = pmodel.p1
  pp = AttractEasyPartnerInterface(p)
  partners.append(pp)
  #Convert peptide
  p = pmodel.p2
  pp = AttractEasyPartnerInterface.empty()
  pp.pdbfile = Data(buildpeptide(p.sequence))
  pp.moleculetype = "Peptide"
  pp.use_rmsd=p.use_rmsd
  pp.rmsd_pdb=p.rmsd_pdb
  pp.ensemble_size = 3
  if pmodel.p1.haddock_restraints:
    pp.haddock_restraints = HaddockRestraintsInterface (
      passivereslist = range(1,len(p.sequence)+1)
     )
  pp = AttractEasyPartnerInterface(pp)
  partners.append(pp)

  newmodel = AttractEasyModel(
   runname=pmodel.runname,
   partners=partners,
   use_grids=pmodel.use_grids,
   rescore_step=False, #does not seem to work well at the present
   use_iattract=pmodel.use_iattract,
   analyze_interface=pmodel.analyze_interface,
   nstruc_analyze_interface=pmodel.nstruc_analyze_interface,
   clustering=pmodel.clustering,
   min_cluster_size=pmodel.min_cluster_size,
   clustering_cutoff=pmodel.clustering_cutoff,
   calc_lrmsd=pmodel.calc_lrmsd,
   calc_irmsd=pmodel.calc_irmsd,
   calc_fnat=pmodel.calc_fnat,
   nr_collect=pmodel.nr_collect,
   np=pmodel.np,
   completion_tool=pmodel.completion_tool,
   use_gpu=pmodel.use_gpu,
  )
  if not pmodel.use_iattract:
      newmodel.max_analysis = 10000

  return newmodel

def buildpeptide(seq):
  import random
  s = random.randint(0,99999)
  while os.path.exists('/tmp/peptide-'+str(s)):
    s = random.randint(0,99999)

  dir = '/tmp/peptide-'+str(s)
  os.system('mkdir '+dir)
  create_PDB(seq,-139,-135,dir+'/model-1.pdb')
  create_PDB(seq,-57,-47,dir+'/model-2.pdb')
  create_PDB(seq,-78,149,dir+'/model-3.pdb')
  for i in range(1,4):
    os.system('echo "MODEL '+str(i)+'\n" >> '+dir+'/peptide.pdb')
    os.system('cat '+dir+'/model-'+str(i)+'.pdb >> '+dir+'/peptide.pdb')
    os.system('echo ENDMDL >> '+dir+'/peptide.pdb')

  data = open(dir+'/peptide.pdb', "r").read()
  #...
  os.system('rm -rf '+dir)
  return data

def create_PDB(seq,phi,psi,output):
  try:
    import Bio.PDB
  except ImportError:
    raise ImportError("You need to have BioPython installed to convert peptide sequences to structures")
  try:
    from PeptideBuilder import Geometry
    import PeptideBuilder
  except ImportError:
    try:
      import Geometry
      import PeptideBuilder
    except ImportError:
      raise ImportError("You need to have the Python module PeptideBuilder installed to convert peptide sequences to structures")

  geo = Geometry.geometry(seq[0])
  struc = PeptideBuilder.initialize_res(geo)
  for aa in seq[1:]:
    struc = PeptideBuilder.add_residue(struc,aa,phi,psi)

  out = Bio.PDB.PDBIO()
  out.set_structure(struc)
  out.save(output)
