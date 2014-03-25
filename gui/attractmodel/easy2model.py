def easy2model(emodel):
  partners = []
  for p in emodel.partners:
    pp = AttractPartnerInterface.empty()
    pp.pdbfile=p.pdbfile
    pp.is_reduced=False
    pp.collect_pdb=p.pdbfile 
    pp.chain="All"
    pp.generate_modes=p.generate_modes
    if p.generate_modes:
      pp.nr_modes=p.nr_modes
    pp.deflex=True
    pp.rmsd_pdb=p.rmsd_pdb
    pp.rmsd_bb=p.rmsd_bb
    if p.ensemble_size > 0:
      pp.ensemble = True
      pp.ensemble_size = p.ensemble_size
      pp.ensemblize = "random"
    pp = AttractPartnerInterface(pp)
    partners.append(pp)
    
  if emodel.use_grids and partners[0].nr_modes > 0:
    partners[0].gridname = "receptorgrid"
    partners[1].gridname = "ligandgrid"
    rgrid = [AttractGrid(gridname="receptorgrid"),AttractGrid(gridname="ligandgrid")]
    iter = [AttractIteration(vmax=1000)]
  
  elif emodel.use_grids:
    partners[0].gridname = "receptorgrid"
    rgrid = [AttractGrid(gridname="receptorgrid")]
    iter = [AttractIteration(vmax=1000)]
    
  else:
    rgrid = []
    iter = [AttractIteration(vmax=50),AttractIteration(rcut=500,vmax=60),AttractIteration(rcut=100,vmax=60),
	    AttractIteration(rcut=50),AttractIteration(rcut=50)]
   
  gravity = 1 if emodel.gravity else 0
  newmodel = AttractModel(runname="easyrun",partners=partners,grids=rgrid,nr_iterations=len(iter),iterations=iter,fix_receptor=True,search="syst",gravity=gravity,calc_lrmsd=emodel.calc_lrmsd,
	       calc_irmsd=emodel.calc_irmsd,calc_fnat=emodel.calc_fnat,nr_collect=emodel.nr_collect,np=emodel.np)
  
  return newmodel