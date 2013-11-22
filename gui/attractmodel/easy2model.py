def easy2model(emodel):
  pass
  
  partners = []
  for p in emodel.partners:
    pp = AttractPartnerInterface(pdbfile=p.pdbfile,chain="All",generate_modes=p.generate_modes,nr_modes=p.nr_modes,
			       deflex=True,rmsd_pdb=p.rmsd_pdb,rmsd_bb=p.rmsd_bb)
    partners.append(pp)
    
  if emodel.use_grids:
    partners[0].gridname = "receptorgrid"
    rgrid = [AttractGrid(gridname="receptorgrid")]
    iter = [AttractIteration(vmax=1000)]
  else:
    rgrid = []
    iter = [AttractIteration(vmax=50),AttractIteration(rcut=500,vmax=60),AttractIteration(rcut=100,vmax=60),
	    AttractIteration(rcut=50),AttractIteration(rcut=50)]
   
  gravity = 1 if emodel.gravity else 0
  newmodel = AttractModel(partners=partners,grids=rgrid,nr_iterations=len(iter),iterations=iter,fix_receptor=True,search="syst",gravity=gravity,calc_lrmsd=emodel.calc_lrmsd,
	       calc_irmsd=emodel.calc_irmsd,calc_fnat=emodel.calc_fnat,nr_collect=emodel.nr_collect,np=emodel.np)
  
  return newmodel