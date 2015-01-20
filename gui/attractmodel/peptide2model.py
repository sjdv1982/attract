peptide2model_version = "Converted from AttractPeptideModel by peptide2model 1.0"
import os
import Geometry, PeptideBuilder, Bio.PDB
import sys
 
def create_PDB(seq,phi,psi,output):
  geo = Geometry.geometry(seq[0])
  struc = PeptideBuilder.initialize_res(geo)
  for aa in seq[1:]:
    struc = PeptideBuilder.add_residue(struc,aa,phi,psi)
    
  out = Bio.PDB.PDBIO()
  out.set_structure(struc)
  out.save(output)

def buildpeptide(seq):
  os.system('mkdir peptide_new')
  create_PDB(seq,-139,-135,'peptide_new/model-1.pdb')
  create_PDB(seq,-57,-47,'peptide_new/model-2.pdb')
  create_PDB(seq,-78,149,'peptide_new/model-3.pdb')
  os.system('python $ATTRACTTOOLS/joinmodel.py peptide_new/model-1.pdb peptide_new/model-2.pdb peptide_new/model-3.pdb > peptide.pdb')
  return 'peptide.pdb'

def peptide2model(pmodel):
  partners = []
  use_flex = []
  for p in pmodel.partners:
    partner_use_flex = False
    pp = AttractPartnerInterface.empty()
    if p.pdbfile is not None:
      pp.pdbfile=p.pdbfile
    else:
      pp.pdbfile = buildpeptide(p.sequence)
      p.ensemble_size = 3
      pp.is_peptide = True
      
    pp.is_reduced=False
    pp.collect_pdb=pp.pdbfile 
    pp.chain="All"
    pp.generate_modes=p.generate_modes
    if p.generate_modes:
      partner_use_flex = True
      pp.nr_modes=p.nr_modes
    pp.deflex=True    
    pp.rmsd_pdb=p.rmsd_pdb
    pp.rmsd_bb=p.rmsd_bb
    if p.ensemble_size > 1:
      partner_use_flex = True
      pp.ensemble = True      
      pp.ensemble_size = p.ensemble_size
      if p.sequence is not None and p.ensemble_size <= 5:
	pp.ensemblize = "all"
      else:
	pp.ensemblize = "random"
    pp = AttractPartnerInterface(pp)
    partners.append(pp)
    use_flex.append(partner_use_flex)
  
  if use_flex[0]:
    receptorgrid = AttractGrid(gridname="receptorgrid", plateau_distance = 10.0, neighbour_distance=12.0)
  else:  
    receptorgrid = AttractGrid(gridname="receptorgrid", plateau_distance = 5.0, neighbour_distance=7.0)

  if use_flex[1]:
    ligandgrid = AttractGrid(gridname="ligandgrid", plateau_distance = 10.0, neighbour_distance=12.0)  
  else:  
    ligandgrid = AttractGrid(gridname="ligandgrid", plateau_distance = 5.0, neighbour_distance=7.0)  
    
  if pmodel.use_grids and partners[0].nr_modes > 0:
    partners[0].gridname = "receptorgrid"
    partners[1].gridname = "ligandgrid"
    rgrid = [receptorgrid,ligandgrid]
    iter = [AttractIteration(vmax=1000)]
  
  elif pmodel.use_grids:
    partners[0].gridname = "receptorgrid"
    rgrid = [receptorgrid]
    iter = [AttractIteration(vmax=1000)]
    
  else:
    rgrid = []
    iter = [
     AttractIteration(vmax=50),
     AttractIteration(rcut=500,vmax=60),
     AttractIteration(rcut=100,vmax=60),
     AttractIteration(rcut=50),
     AttractIteration(rcut=50),
    ]
   
  gravity = 2 if pmodel.gravity else 0
  newmodel = AttractModel (
   runname=pmodel.runname,
   partners=partners,
   grids=rgrid,
   nr_iterations=len(iter),
   iterations=iter,
   fix_receptor=True,
   search="random",
   structures=100000,
   gravity=gravity,
   calc_lrmsd=pmodel.calc_lrmsd,
   calc_irmsd=pmodel.calc_irmsd,
   calc_fnat=pmodel.calc_fnat,
   nr_collect=pmodel.nr_collect,
   np=pmodel.np,
   deredundant_ignorens = True,
   annotation = peptide2model_version,
  )  
  if pmodel.use_iattract:
    iattract = IAttractParameters(
     nstruc = pmodel.nr_collect,
     icut = 5.0
    ) 
    newmodel.iattract = iattract
    newmodel.demode = True ##TODO: maybe we want this in all cases..???
    newmodel.validate()
  return newmodel