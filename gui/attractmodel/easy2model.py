easy2model_version = "Converted from AttractEasyModel by easy2model 1.0"

def easy2model(emodel):
  partners = []
  use_flex = []
  for p in emodel.partners:
    partner_use_flex = False
    pp = AttractPartnerInterface.empty()
    pp.pdbfile=p.pdbfile
    pp.is_reduced=False
    pp.collect_pdb=p.pdbfile 
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
    
  if emodel.use_grids and partners[0].nr_modes > 0:
    partners[0].gridname = "receptorgrid"
    partners[1].gridname = "ligandgrid"
    rgrid = [receptorgrid,ligandgrid]
    iter = [AttractIteration(vmax=1000)]
  
  elif emodel.use_grids:
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
   
  gravity = 2 if emodel.gravity else 0
  newmodel = AttractModel (
   runname=emodel.runname,
   partners=partners,
   grids=rgrid,
   nr_iterations=len(iter),
   iterations=iter,
   fix_receptor=True,
   search="syst",
   gravity=gravity,
   calc_lrmsd=emodel.calc_lrmsd,
   calc_irmsd=emodel.calc_irmsd,
   calc_fnat=emodel.calc_fnat,
   nr_collect=emodel.nr_collect,
   np=emodel.np,
   deredundant_ignorens = False,
   annotation = easy2model_version,
  )
  return newmodel