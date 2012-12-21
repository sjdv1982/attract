import os

#datamodel: see TODOs
#TODO: add check between p.modes_file and p.nr_modes
#TODO: add check that either all interactions are grid-accelerated, or none are
#  (and an option to disable this check; for this, adapt --rcut iteration as well)
#TODO: a script to add energies back to deredundant output

def generate(m):
  if m.forcefield != "ATTRACT": 
    reduced_all = all((p.is_reduced for p in m.partners))
    if not reduced_all:
      #TODO?
      raise ValueError(
"""
When not using the ATTRACT forcefield, the 'reduce' program may easily give errors
Therefore, you can only define docking partners that are already in the reduced form
"""
      ) 
  
  if m.forcefield == "ATTRACT":
    ffpar = "$ATTRACTDIR/../parmw.par"
  elif m.forcefield == "OPLSX":  
    ffpar = "$ATTRACTDIR/../allatom/allatom.par"

  ret = ""
  cleanupfiles = []
  if (m.header is not None):
    ret = m.header + "\n\n"
  ret += "if [ 1 -eq 1 ]; then ### move and change to disable parts of the protocol\n"
  ret += "$ATTRACTDIR/shm-clean\n\n"

  modes_any = any((p.modes_file for p in m.partners))
  aa_modes_any = any((p.aa_modes_file for p in m.partners))
  if modes_any:
    ret += """
echo '**************************************************************'
echo 'Assemble modes file...'
echo '**************************************************************'
cat /dev/null > hm-all.dat
""" 
    if aa_modes_any:
      hm_all_aa = "hm-all-aa.dat"
      ret += "cat /dev/null > hm-all-aa.dat\n"
    else: 
      hm_all_aa = "hm-all.dat"
    ret += "\n"  
    for p in m.partners:
      if p.generate_modes: 
        raise ValueError("Generating modes must currently be done manually") #TODO?
      if p.moleculetype != "Protein": raise Exception #TODO: not yet implemented
      if p.nr_modes is None or p.nr_modes == 0 or  p.modes_file is None:
        ret += "echo 0 >> hm-all.dat\n"
        if aa_modes_any:
          ret += "echo 0 >> hm-all-aa.dat\n"
      elif p.nr_modes == 10:
        ret += "cat %s >> hm-all.dat\n" % p.modes_file.name
        mf = p.modes_file.name
        if p.aa_modes_file is not None: mf = p.aa_modes_file.name
        if aa_modes_any:
          ret += "cat %s >> hm-all-aa.dat\n" % mf
      else:
        ret += "#select first %d modes...\n" % (p.nr_modes)
        ret += "awk 'NF == 2 && $1 == %d{exit}{print $0}' %s >> hm-all.dat\n" % \
         (p.nr_modes+1, p.modes_file.name)
        if aa_modes_any:
          mf = p.modes_file.name
          if p.aa_modes_file is not None: mf = p.aa_modes_file.name
          ret += "awk 'NF == 2 && $1 == %d{exit}{print $0}' %s >> hm-all-aa.dat\n" % \
           (p.nr_modes+1, mf)         
    ret += "\n"    
  filenames = []
  reduce_any = False
  reduced = set()
  for pnr,p in enumerate(m.partners):
    if p.pdb.mode == "download":
      raise Exception("TODO: not yet implemented")
    if p.is_reduced == False:
      if reduce_any == False:
        ret += """
echo '**************************************************************'
echo 'Reduce partner PDBs...'
echo '**************************************************************'
"""      
        reduce_any = True
      pdbname = p.pdb.pdbfile.name
      pdbname2 = os.path.split(pdbname)[1]
      pdbname_reduced = pdbname2[:-4] + "r.pdb"
      if pdbname_reduced not in reduced:
        if pdbname2 != pdbname:
          ret += "cat %s > %s\n" % (pdbname, pdbname2)
        ret += "$ATTRACTDIR/reduce %s > /dev/null\n" % pdbname2
        reduced.add(pdbname_reduced)
      filenames.append(pdbname_reduced)
    else:  
      filenames.append(p.pdb.pdbfile.name) 
    #TODO?: generate normal modes  
  if reduce_any: ret += "\n"
  
  if len(m.partners) == 2:
    partnerfiles = " ".join(filenames)
  else:
    partnerfiles = "partners.pdb"
    ret += """
echo '**************************************************************'
echo 'Concatenate partner PDBs...'
echo '**************************************************************'
echo %d > partners.pdb
""" % len(m.partners)
    for f in filenames:
      ret += "grep ATOM %s >> partners.pdb\n" % f
      ret += "echo TER >> partners.pdb\n"
          
  ret += """
#name of the run
name=%s 
""" % m.runname
  params = "\"" + ffpar + " " + partnerfiles
  scoreparams = params + " --score --fix-receptor"
  gridparams = ""
  if m.fix_receptor: params += " --fix-receptor" 
  if modes_any: 
    ps = " --modes hm-all.dat"
    params += ps
    scoreparams += ps
  gridfiles = {}
  ret_shm = ""
  for g in m.grids:
    if g.gridfile is not None: 
      v = g.gridfile.name
      if m.np > 1:
        v2 = v +"header"      
        cleanupfiles.append(v2)
        tq = "-torque" if g.torque else ""
        ret_shm = "\n#Load %s grid into memory\n" % (g.gridname)
        ret_shm += "$ATTRACTDIR/shm-grid%s %s %s\n" % (tq, v, v2)
        v = v2
    else:  
      gheader = ""
      if m.np > 1: gheader = "header"
      v = g.gridname.strip() + ".grid" + gheader
    gridfiles[g.gridname.strip()] = v
  grid_used = {}
  ens_any = False
  for pnr,p in enumerate(m.partners):
    if p.gridname is not None:
      v = gridfiles[p.gridname.strip()]
      if v in grid_used: 
        v = grid_used[v]
      else:
        grid_used[v] = pnr+1
      gridparams += " --grid %d %s" % (pnr+1, str(v))
    if p.ensemble_list is not None  :
      ps = " --ens %d %s" % (pnr+1, p.ensemble_list.name)
      params += ps
      scoreparams += ps
      ens_any = True
  if m.ghost:
    params += " --ghost"
  if m.gravity:
    params += " --gravity %d" % m.gravity
  if m.rstk is not None and m.rstk != 0.2:
    params += " --rstk %s" % str(m.rstk)
  if m.dielec == "cdie": 
    ps = " --cdie"
    params += ps  
    scoreparams += ps  
  if m.epsilon != 15: 
    ps = " --epsilon %s" % (str(m.epsilon))
    params += ps
    scoreparams += ps  
    
  for sym in m.symmetries:
    symcode = len(sym.partners)
    if sym.symmetry == "Dx": symcode = -4
    partners = " ".join([str(v) for v in sym.partners])
    params += " --sym %d %s" % (symcode, partners)
  if m.cryoem_data:
    params += " --em %s" % m.cryoem_data.name
  params += "\""
  scoreparams += "\""
  ret += """
#docking parameters
params=%s
scoreparams=%s
"""  % (params, scoreparams)
  if len(gridparams):
    ret += """
#grid parameters
gridparams="%s"
"""  % gridparams
  
  if m.np > 1:
    if m.jobsize in (0, None):
      parals = "--np %d --chunks %d" % (m.np, m.np)
    else:
      parals = "--np %d --jobsize %d" % (m.np, m.jobsize)
    ret += """
#parallelization parameters
parals="%s"
"""  % (parals)
  
  ret += ret_shm
  if m.search == "syst" or m.search == "custom":
    if m.search == "syst" or m.start_structures_file is None:
      ret += """
echo '**************************************************************'
echo 'generate starting structures...'
echo '**************************************************************'
"""
      rotfil = "$ATTRACTDIR/../rotation.dat"
      if m.rotations_file is not None:
        rotfil = m.rotations_file.name
      if rotfil != "rotation.dat":
        ret += "cat %s > rotation.dat\n" % rotfil
      if m.translations_file is not None:
        transfil = m.translations_file.name
        if transfil != "translate.dat":
          ret += "cat %s > translate.dat\n" % transfil
      else:
        ret += "$ATTRACTDIR/translate %s %s > translate.dat\n" % \
         (filenames[0], filenames[1])
      ret += "$ATTRACTDIR/systsearch > systsearch.dat\n"
      ret += "start=systsearch.dat\n"
      start = "systsearch.dat"
    else:
      ret += """
#starting structures
start=%s 
""" % m.start_structures_file.name
      start = m.start_structures_file.name
  elif m.search == "random":
    ret += "python $ATTRACTTOOLS/randsearch.py %d %d > randsearch.dat\n" % \
     (len(m.partners), m.structures)
    ret += "start=randsearch.dat\n"    
    start = randsearch.dat
  else:
    raise Exception("Unknown value")
  ret += "\n"  
  
  inp = "$start"
  if any((p.ensemblize for p in m.partners if p.ensemblize not in (None,"custom"))):
    start0 = start
    ret += """
echo '**************************************************************'
echo 'ensemble search:' 
echo ' add ensemble conformations to the starting structures'
echo '**************************************************************'
"""
    for pnr, p in enumerate(m.partners):
      if p.ensemblize in (None, "custom"): continue
      if p.ensemblize == "all":
        ret += """
echo '**************************************************************'
echo 'multiply the number of structures by the ensemble for ligand %d'
echo '**************************************************************'
""" % (pnr+1)
      elif p.ensemblize == "random":
        ret += """
echo '**************************************************************'
echo 'random ensemble conformation in ligand %d for each starting structure'
echo '**************************************************************'
""" % (pnr+1)
      else: 
        raise ValueError(p.ensemblize)
      if start == start0 and start not in ("systsearch.dat", "randsearch.dat"):
        start2 = "start-ens%d.dat" % (pnr+1)
      else:
        start2 = os.path.splitext(start)[0] + "-ens%d.dat" % (pnr+1)
      ret += "python $ATTRACTTOOLS/ensemblize.py %s %d %d %s > %s\n" % \
       (start, p.ensemble_size, pnr+1, p.ensemblize, start2)
      start = start2
      inp = start 
    
  for g in m.grids:
    if g.gridfile is not None: continue
    gridfile = gridfiles[g.gridname.strip()]
    if gridfile not in grid_used: continue
    ret += """
echo '**************************************************************'
echo 'calculate %s grid'
echo '**************************************************************'
""" % g.gridname
    no = ""
    if not g.mask_interior: no = "no_"
    partner = grid_used[gridfile]-1
    f = filenames[partner]
    f0 = os.path.splitext(f)[0]
    vol = "%s-%sinterior.vol" % (f0, no)
    ret += "$ATTRACTDIR/calc_%sinterior %s %s\n" % (no, f, vol)    
    tomp = ""
    if g.torque: tomp += "-torque"
    if g.omp: tomp += "-omp"
    tail = ""
    if m.np > 1: 
      tail += " --shm"
    if m.forcefield == "OPLSX":  
      tail += " --alphabet $ATTRACTDIR/../allatom/oplsx.trans"
    if m.dielec == "cdie": tail += " --cdie"
    if m.epsilon != 15: tail += " --epsilon %s" % (str(m.epsilon))    
    if g.calc_potentials == False: tail += " -calc-potentials=0"
    if g.omp:
      ret += "for i in `seq 1 10`; do\n\n"      
    ret += "$ATTRACTDIR/make-grid%s %s %s %s %s %s %s %s\n\n" % \
     (tomp, f, vol, ffpar, g.plateau_distance, g.neighbour_distance, gridfile, tail)
    if g.omp:
      ret += """
if [ $? = 0 ]; then
  break
fi
echo 'Retry grid calculation...'
done
"""
  
  ret += """
echo '**************************************************************'
echo 'Docking'
echo '**************************************************************'
"""
  ordinals = ["1st", "2nd", "3rd",] + ["%dth" % n for n in range(4,51)]
  iterations = []
  for n in range(m.nr_iterations):  
    if m.iterations is None or len(m.iterations) <= n:
      iterations.append([None, None, False, False])
    else:
      it = m.iterations[n]
      newit = [it.rcut, it.vmax, it.traj, it.mc]
      if it.mc:
        mcpar = (it.mctemp, it.mcscalerot, it.mcscalecenter, it.mcscalemode, it.mcensprob)
        newit.append(mcpar)
      iterations.append(newit)
  for i,it in enumerate(iterations):
    ret += """
echo '**************************************************************'
echo '%s minimization'
echo '**************************************************************'
""" % ordinals[i]
    itparams = ""
    rcut, vmax, traj, mc = it[:4]
    if rcut is not None and len(grid_used) == 0: itparams += " --rcut %s" % str(rcut)
    if vmax is not None: itparams += " --vmax %s" % str(vmax)
    if traj: itparams += " --traj"
    if mc: 
      itparams += " --mc"
      mctemp, mcscalerot, mcscalecenter, mcscalemode, mscensprob = it[4]
      if mctemp is not None: itparams += " --mctemp %s" % str(mctemp)
      if mcscalerot is not None: itparams += " --mcscalerot %s" % str(mcscalerot)
      if mcscalecenter is not None: itparams += " --mcscalecenter %s" % str(mcscalecenter)
      if mcensprob is not None: itparams += " --mcensprob %s" % str(mcensprob)
    if m.np > 1:
      attract = "python $ATTRACTDIR/../protocols/attract.py"
      tail = "$parals --output"  
    else:
      attract = "$ATTRACTDIR/attract"
      tail = ">"  
    if i == len(iterations) - 1:
      outp = "out_$name.dat"
    else:  
      outp = "stage%d_$name.dat" % (i+1)  
    gridpar = ""
    if len(gridparams): gridpar = " $gridparams" 
    ret += "%s %s $params%s %s %s %s\n" % (attract, inp, gridpar, itparams, tail, outp)
    inp = outp       
  ret += "\n"  

  result = outp  
  if m.rescoring:
    ret += """
echo '**************************************************************'
echo 'Final rescoring'
echo '**************************************************************'
"""
    if m.np > 1:
      attract = "python $ATTRACTDIR/../protocols/attract.py"
      tail = "$parals --output"  
    else:
      attract = "$ATTRACTDIR/attract"
      tail = ">"  
    ret += "%s %s $scoreparams --rcut %s %s out_$name.score\n" \
     % (attract, result, str(m.rcut_rescoring), tail)
    ret += """     
echo '**************************************************************'
echo 'Merge the scores with the structures'
echo '**************************************************************'
"""
    ret += "python $ATTRACTTOOLS/fill-energies.py %s out_$name.score > out_$name-scored.dat\n" % (result)
    ret += "\n"  
    result = "out_$name-scored.dat"

  if m.sort:
    ret += """
echo '**************************************************************'
echo 'Sort structures'
echo '**************************************************************'
"""
    ret += "python $ATTRACTTOOLS/sort.py %s > out_$name-sorted.dat\n" % result
    ret += "\n"  
    result = "out_$name-sorted.dat"

  if m.deredundant:
    outp = os.path.splitext(result)[0] +"-dr.dat"
    ret += """
echo '**************************************************************'
echo 'Remove redundant structures'
echo '**************************************************************'
""" 
    par_flex = ""
    if ens_any:
      par_flex += " --ens"
      for p in m.partners:
        ensemble_size = p.ensemble_size
        if ensemble_size is None: ensemble_size = 0
        par_flex += " %d" % ensemble_size
      if m.deredundant_ignorens: par_flex += " --ignorens"
    if modes_any:
      par_flex += " --modes"
      for p in m.partners:
        nr_modes = p.nr_modes
        if nr_modes is None: nr_modes = 0
        par_flex += " %d" % nr_modes
    
    ret += "$ATTRACTDIR/deredundant %s %d%s > %s\n" % (result, len(m.partners), par_flex, outp)
    ret += "\n"  
    result = outp
  result0 = result
  if m.calc_lrmsd or m.calc_irmsd:
    deflex_any = any((p.deflex for p in m.partners))
    if deflex_any:
      ret += """
echo '**************************************************************'
echo 'Remove flexibility for RMSD calculations'
echo '**************************************************************'
tmpf=`mktemp`
tmpf2=`mktemp`

""" 
      outp, outp2 = "$tmpf", "$tmpf2"
      for pnr,p in enumerate(m.partners):
        if p.deflex:
          if p.nr_modes:           
            ret += "python $ATTRACTTOOLS/demode.py %s %d > %s\n" % \
             (result, pnr+1,outp)
          elif p.ensemble_size:
            ret += "python $ATTRACTTOOLS/de-ensemblize.py %s %d > %s\n" % \
             (result, pnr+1,outp)
          else:
            continue
          result = outp
          outp, outp2 = outp2, outp
  if m.calc_lrmsd or m.calc_irmsd:    
    any_bb = any((p.rmsd_bb for p in m.partners[1:]))
    if any_bb:
      ret += """
echo '**************************************************************'
echo 'Select backbone atoms'
echo '**************************************************************'
""" 
    rmsd_filenames = []
    for pnr in range(len(m.partners)):
      filename = filenames[pnr]
      p = m.partners[pnr]
      if p.rmsd_bb:
        filename2 = os.path.splitext(filename)[0] + "-bb.pdb"
        ret += "$ATTRACTTOOLS/backbone %s > %s\n" % (filename, filename2)
        filename = filename2                
      rmsd_filenames.append(filename)
    if m.calc_irmsd or m.calc_fnat:
      irmsd_filenames = rmsd_filenames[:]
      irmsd_refenames = rmsd_filenames[:]
      for pnr in range(len(m.partners)):
	filename = filenames[pnr]
	p = m.partners[pnr]
	if p.rmsd_pdb is not None:
          filename = p.rmsd_pdb.name
          if p.rmsd_bb:
            filename2 = os.path.splitext(filename)[0] + "-bb.pdb"
            ret += "$ATTRACTTOOLS/backbone %s > %s\n" % (filename, filename2)
            filename = filename2                
	elif p.rmsd_bb:
          filename2 = os.path.splitext(filename)[0] + "-bb.pdb"
          filename = filename2          
	irmsd_refenames[pnr] = filename
    if m.calc_lrmsd:
      lrmsd_filenames = rmsd_filenames[1:]
      if m.calc_irmsd or m.calc_fnat:
        lrmsd_refenames = irmsd_refenames[1:]
      else:
        lrmsd_refenames = rmsd_refenames[1:]
	for pnr in range(1, len(m.partners)):
	  filename = filenames[pnr]
	  p = m.partners[pnr]
	  if p.rmsd_pdb is not None:
            filename = p.rmsd_pdb.name
            if p.rmsd_bb:
              filename2 = os.path.splitext(filename)[0] + "-bb.pdb"
              ret += "$ATTRACTTOOLS/backbone %s > %s\n" % (filename, filename2)
              filename = filename2                
	  elif p.rmsd_bb:
            filename2 = os.path.splitext(filename)[0] + "-bb.pdb"
            filename = filename2          
	  lrmsd_refenames[pnr-1] = filename
    if any_bb:  
      ret += "\n"
      bb_str = "backbone "
    else:
      bb_str = ""
  if m.calc_lrmsd or m.calc_irmsd or m.calc_fnat: 
    flexpar = ""
    if modes_any and not deflex_any: 
      if aa_modes_any:
        flexpar = " --modes hm-all-aa.dat"
      else:
        flexpar = " --modes hm-all.dat"
    for pnr,p in enumerate(m.partners):
      if p.deflex == False and p.ensemble_list is not None:
        flexpar += " --ens %d %s" % (pnr+1,p.ensemble_list.name)
  if m.calc_lrmsd:      
    ret += """
echo '**************************************************************'
echo 'calculate %sligand RMSD'
echo '**************************************************************'
""" % bb_str      
	
    fixresult = None 
    if m.fix_receptor == False: 
      fixresult = result + "-fixre"
      ret += "$ATTRACTDIR/fix_receptor %s %d%s > %s\n" % (result, len(m.partners), flexpar, result2)
      result = fixresult
	
    lrmsd_allfilenames = []
    for f1, f2 in zip(lrmsd_filenames, lrmsd_refenames):
      lrmsd_allfilenames.append(f1)
      lrmsd_allfilenames.append(f2)
    lrmsd_allfilenames = " ".join(lrmsd_allfilenames)
    lrmsdresult = os.path.splitext(result0)[0] + ".lrmsd"
    ret += "$ATTRACTDIR/lrmsd %s %s%s > %s\n" % (result, lrmsd_allfilenames, flexpar, lrmsdresult)
    ret += "\n"

  if m.calc_irmsd:      
    ret += """
echo '**************************************************************'
echo 'calculate %sinterface RMSD'
echo '**************************************************************'
""" % bb_str      
		
    irmsd_allfilenames = []
    for f1, f2 in zip(irmsd_filenames, irmsd_refenames):
      irmsd_allfilenames.append(f1)
      irmsd_allfilenames.append(f2)
    irmsd_allfilenames = " ".join(irmsd_allfilenames)
    irmsdresult = os.path.splitext(result0)[0] + ".irmsd"
    ret += "python $ATTRACTDIR/irmsd.py %s %s%s > %s\n" % (result, irmsd_allfilenames, flexpar, irmsdresult)
    ret += "\n"

  if m.calc_lrmsd or m.calc_irmsd:
    if deflex_any:
      ret += "rm -f $tmpf $tmpf2\n"
      result = result0
    if m.fix_receptor == False:
      ret += "rm -f %s\n" % fixresult

  if m.calc_fnat:      
    ret += """
echo '**************************************************************'
echo 'calculate fraction of native contacts'
echo '**************************************************************'
"""		
    fnat_allfilenames = []
    for f1, f2 in zip(irmsd_filenames, irmsd_refenames):
      fnat_allfilenames.append(f1)
      fnat_allfilenames.append(f2)
    fnat_allfilenames = " ".join(fnat_allfilenames)
    fnatresult = os.path.splitext(result0)[0] + ".fnat"
    ret += "python $ATTRACTDIR/fnat.py %s 5 %s%s > %s\n" % (result, fnat_allfilenames, flexpar, fnatresult)
    ret += "\n"

  if m.collect:
    collect_filenames = filenames
    nr = m.nr_collect
    for pnr in range(len(m.partners)):
      p = m.partners[pnr]
      if p.collect_pdb is not None:
        collect_filenames[pnr] = p.collect_pdb.name
    ret += """
echo '**************************************************************'
echo 'collect top %d structures:'
echo '**************************************************************'
""" % nr
    ret += "$ATTRACTTOOLS/top %s %d > out_$name-top%d.dat\n" % (result, nr, nr)
    collect_filenames = " ".join(collect_filenames)
    flexpar_aa = ""
    if modes_any: 
      if aa_modes_any:
        flexpar_aa = " --modes hm-all-aa.dat"
      else:
        flexpar_aa = " --modes hm-all.dat"    
    for pnr,p in enumerate(m.partners):
      if p.collect_ensemble_list is not None:
        flexpar_aa += " --ens %d %s" % (pnr+1,p.collect_ensemble_list.name)    
      elif p.ensemble_list is not None:
        flexpar_aa += " --ens %d %s" % (pnr+1,p.ensemble_list.name)    
    ret += "$ATTRACTDIR/collect out_$name-top%d.dat %s%s > out_$name-top%d.pdb\n" % \
     (nr, collect_filenames, flexpar_aa, nr)
    ret += "\n"
  if len(cleanupfiles):  
    ret += "\nrm -f " + " ".join(cleanupfiles) + "\n"
  if (m.footer is not None):
    ret += m.footer + "\n\n"
    
  ret += "fi ### move to disable parts of the protocol\n"
  ret = ret.replace("\n\n\n","\n\n")
  ret = ret.replace("\n\n\n","\n\n")
  return ret

