easy2model_version = "Converted from AttractEasyModel by easy2model 1.1"

def easy2model(emodel):
  partners = []
  use_flex = []
  use_haddock = False
  has_peptide = False
  has_na = False
  for i,p in enumerate(emodel.partners):
    partner_use_flex = False
    pp = AttractPartnerInterface.empty()
    pp.pdbfile=p.pdbfile
    pp.has_hydrogens=p.has_hydrogens
    pp.is_reduced=False
    #pp.collect_pdb=p.pdbfile 
    pp.chain="All"
    pp.haddock_restraints = p.haddock_restraints
    if p.haddock_restraints: use_haddock = True
    if p.moleculetype == "Peptide":
      pp.moleculetype = "Protein"
      pp.charged_termini = True
      has_peptide = True
    else:
      pp.moleculetype=p.moleculetype
      
    pp.generate_modes=p.generate_modes
    if p.generate_modes:
      partner_use_flex = True
      pp.nr_modes=p.nr_modes
    if p.moleculetype in ("DNA", "RNA"):  
      pp.aacontact_modes = True
      has_na = True
    if emodel.use_iattract:
      pp.deflex=False
    else:
      pp.deflex=True    
    
    pp.rmsd_pdb=p.rmsd_pdb
    if p.ensemble_size > 1:
      partner_use_flex = True
      pp.ensemble = True      
      pp.ensemble_size = p.ensemble_size
      if p.moleculetype == "Peptide":
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
  
  rgrid = []
  if emodel.use_grids:
    if partners[0].nr_modes > 0:
      partners[0].gridname = "receptorgrid"
      partners[1].gridname = "ligandgrid"
      rgrid = [receptorgrid,ligandgrid]
    else:
      partners[0].gridname = "receptorgrid"
      rgrid = [receptorgrid]      
  
  if use_haddock:
    iter = [
     AttractIteration(vmax=50, prep=True),    
     AttractIteration(vmax=1000),
     AttractIteration(vmax=1000, restweight=0.01),
    ]    
  elif emodel.use_grids:    
    iter = [AttractIteration(vmax=1000)]    
  else:    
    iter = [
     AttractIteration(vmax=50),
     AttractIteration(rcut=500,vmax=60),
     AttractIteration(rcut=100,vmax=60),
     AttractIteration(rcut=50),
     AttractIteration(rcut=50),
    ]
   
  gravity = 2
  if use_haddock: gravity = 0
  if not has_peptide and not has_na: gravity = 0
  rmsd_atoms = "backbone"
  if has_na: 
    rmsd_atoms = "all"
  newmodel = AttractModel (
   runname=emodel.runname,
   partners=partners,
   grids=rgrid,
   nr_iterations=len(iter),
   iterations=iter,
   fix_receptor=True,
   search="syst",
   gravity=gravity,
   forcefield=emodel.forcefield,
   calc_lrmsd=emodel.calc_lrmsd,
   calc_irmsd=emodel.calc_irmsd,
   rmsd_atoms=rmsd_atoms,
   calc_fnat=emodel.calc_fnat,
   nr_collect=emodel.nr_collect,
   np=emodel.np,
   deredundant_ignorens = False,
   demode=True,
   completion_tool=emodel.completion_tool,
   annotation = easy2model_version,   
  )  
  if has_peptide or has_na or use_haddock:
    newmodel.search = "random"
    if has_na:
      newmodel.structures= 200000
    else:
      newmodel.structures= 100000
    
  if emodel.use_iattract:
    iattract = IAttractParameters(
     nstruc = emodel.nr_collect
    )
    if has_peptide:
      iattract.icut= 5.0      
    newmodel.iattract = iattract
    newmodel.validate()
  return newmodel