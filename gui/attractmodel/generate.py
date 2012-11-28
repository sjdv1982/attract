import os

def generate(m):
  if m.forcefield != "ATTRACT": raise Exception #TODO, not yet implemented
  if m.calc_irmsd: raise Exception #TODO, not yet implemented
  ret = ""
  #ret += "#!/bin/sh\n"
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
      if p.generate_modes: raise Exception #TODO: not yet implemented
      if p.moleculetype != "Protein": raise Exception #TODO: not yet implemented
      if p.nr_modes is None or p.nr_modes == 0:
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
        ret += "$ATTRACTDIR/reduce %s >& /dev/null\n" % pdbname2
        reduced.add(pdbname_reduced)
      filenames.append(pdbname_reduced)
    else:  
      filenames.append(p.pdb.pdbfile.name) 
    #TODO: generate normal modes  
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
  params = "\"$ATTRACTDIR/../parmw.par " + partnerfiles
  gridparams = ""
  if m.fix_receptor: params += " --fix-receptor" 
  if modes_any: params += " --modes hm-all.dat"
  gridfiles = {}
  ret_shm = ""
  for g in m.grids:
    if g.gridfile is not None: 
      v = g.gridfile.name
      if m.np > 1:
        v2 = v +"header"      
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
      params += " --ens %d %s" % (pnr+1, p.ensemble_list.name)
      ens_any = True
  if m.ghost:
    params += " --ghost"
  if m.gravity:
    params += " --gravity %d" % m.gravity
  if m.rstk is not None and m.rstk != 0.2:
    params += " --rstk %s" % str(m.rstk)
    
  for sym in m.symmetries:
    symcode = len(sym.partners)
    if sym.symmetry == "Dx": symcode = -4
    partners = " ".join([str(v) for v in sym.partners])
    params += " --sym %d %s" % (symcode, partners)
  if m.cryoem_data:
    params += " --em %s" % m.cryoem_data.name
  params += "\""
  ret += """
#docking parameters
params=%s
"""  % params
  if len(gridparams):
    ret += """
#grid parameters
gridparams="%s"
"""  % gridparams
  
  if m.np > 1:
    ret += """
#parallelization parameters
parals="--jobsize %d --np %d"
"""  % (m.jobsize, m.np)
  
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
    if m.dielec == "cdie": tail += " --cdie"
    if m.epsilon != 15: tail += " --epsilon %s" % (str(g.epsilon))
    #TODO: non-ATTRACT forcefield
    if g.calc_potentials == False: tail += " -calc-potentials=0"
    if g.omp:
      ret += "for i in `seq 1 10`; do\n\n"      
    ret += "$ATTRACTDIR/make-grid%s %s %s $ATTRACTDIR/../parmw.par %s %s %s %s\n\n" % \
     (tomp, f, vol, g.plateau_distance, g.neighbour_distance, gridfile, tail)
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
    if rcut is not None and len(grid_used) > 0: itparams += " --rcut %s" % str(rcut)
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
    ret += "$ATTRACTDIR/attract %s $params --rcut %s --score > out_$name.score\n" \
     % (result, str(m.rcut_rescoring))
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
    par_ens = ""
    if ens_any:
      par_ens = " --ens"
      for p in m.partners:
        ensemble_size = p.ensemble_size
        if ensemble_size is None: ensemble_size = 1
        par_ens += " %d" % ensemble_size
      if m.deredundant_ignorens: par_ens += " --ignorens"
    ret += "python $ATTRACTTOOLS/deredundant.py %s%s > %s\n" % (result, par_ens, outp)
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
  if m.calc_lrmsd:
    if m.fix_receptor == False: raise Exception #TODO
    any_ca = any((p.lrmsd_ca for p in m.partners[1:]))
    if any_ca:
      ret += """
echo '**************************************************************'
echo 'Select CA atoms'
echo '**************************************************************'
""" 
    lrmsd_refenames = filenames[1:]
    for pnr in range(1, len(m.partners)):
      filename = filenames[pnr]
      p = m.partners[pnr]
      if p.rmsd_pdb is not None:
        filename = p.rmsd_pdb.name
        if p.lrmsd_ca:
          filename2 = os.path.splitext(filename)[0] + "ca.pdb"
          ret += "grep ' CA ' %s > %s\n" % (filename, filename2)
          filename = filename2                
      elif p.lrmsd_ca:
        filename2 = os.path.splitext(filename)[0] + "ca.pdb"
        filename = filename2          
      lrmsd_refenames[pnr-1] = filename
  
    lrmsd_filenames = filenames[1:]
    for pnr in range(1, len(m.partners)):
      filename = filenames[pnr]
      p = m.partners[pnr]
      if p.lrmsd_ca:
        filename2 = os.path.splitext(filename)[0] + "ca.pdb"
        ret += "grep ' CA ' %s > %s\n" % (filename, filename2)
        filename = filename2
      lrmsd_filenames[pnr-1] = filename
    if any_ca:  
      ret += "\n"
      ca_str = "CA "
    else:
      ca_str = ""
    ret += """
echo '**************************************************************'
echo 'calculate %sligand RMSD'
echo '**************************************************************'
""" % ca_str   
    lrmsd_allfilenames = []
    for f1, f2 in zip(lrmsd_filenames, lrmsd_refenames):
      lrmsd_allfilenames.append(f1)
      lrmsd_allfilenames.append(f2)
    lrmsd_allfilenames = " ".join(lrmsd_allfilenames)
    lrmsdresult = os.path.splitext(result0)[0] + ".lrmsd"
    par = ""
    if modes_any and not deflex_any: 
      if aa_modes_any:
        par = " --modes hm-all-aa.dat"
      else:
        par = " --modes hm-all.dat"
    for pnr,p in enumerate(m.partners):
      if p.deflex == False and p.ensemble_list is not None:
        par += " --ens %d %s" % (pnr+1,p.ensemble_list.name)
    ret += "$ATTRACTDIR/lrmsd %s %s%s > %s\n" % (result, lrmsd_allfilenames, par, lrmsdresult)
    ret += "\n"

  if m.calc_lrmsd or m.calc_irmsd:
    if deflex_any:
      ret += "rm -f $tmpf $tmpf2\n"
      result = result0

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
    par = ""
    if modes_any: 
      if aa_modes_any:
        par = " --modes hm-all-aa.dat"
      else:
        par = " --modes hm-all.dat"    
    for pnr,p in enumerate(m.partners):
      if p.collect_ensemble_list is not None:
        par += " --ens %d %s" % (pnr+1,p.collect_ensemble_list.name)    
      elif p.ensemble_list is not None:
        par += " --ens %d %s" % (pnr+1,p.ensemble_list.name)    
    ret += "$ATTRACTDIR/collect out_$name-top%d.dat %s%s > out_$name-top%d.pdb\n" % \
     (nr, collect_filenames, par, nr)
    ret += "\n"
  
    
  ret += "fi ### move to disable parts of the protocol\n"
  ret = ret.replace("\n\n\n","\n\n")
  ret = ret.replace("\n\n\n","\n\n")
  return ret

#raise Exception("Not implemented")
"""
def generate(m):
  ret = ""
  if (m.header is not None):
    ret = m.header + "\n\n"
    
  ret += "$ATTRACTDIR/bin/shm-clean\n\n"
        
  npart = len(m.partners)  
  
  if m.search == "Random":
    ret += "python $ATTRACTDIR/tools/randsearch.py %d %d > start.dat\n\n" % (npart, m.structures)
  else:
    raise Exception("TODO")
    
  if m.use_grids:
    ngrid = npart
    if npart == 2 and m.fix_receptor:
      ngrid = 1
    elif m.grid_torque:
      ngrid -= 1  
    if m.grid_mask_interior:
      calc_interior = "$ATTRACTDIR/bin/calc_interior"
    else:
      calc_interior = "$ATTRACTDIR/bin/calc_no_interior"
    for pnr in range(ngrid):
      part = "partner-%d" % (pnr+1)
      ret += calc_interior + " %s.pdb %s-interior.vol\n" % (part, part)
    ret += "\n"

    if m.forcefield == "ATTRACT":
      for pnr in range(npart):    
        part = "partner-%d" % (pnr+1)
        s = "$ATTRACTDIR/bin/reduce %s.pdb > /dev/null" % part
        ret += s + "\n"
      if npart == 2:
        pdbfiles = "partner-1r.pdb partner-2r.pdb"
      else:
        ret += "echo %d > partners.pdb\n" % npart
        for pnr in range(npart):
          ret += "grep ATOM partner-%dr.pdb >> partners.pdb\n" % (pnr+1)
          ret += "echo TER >> partners.pdb\n" 
        ret += "echo END >> partners.pdb\n"
        pdbfiles = "partners.pdb"
      ret += "\n"  
          
    else:
      raise Exception("TODO")
       
    make_grid = "$ATTRACTDIR/bin/make-grid"
    shm_grid = "$ATTRACTDIR/bin/shm-grid" 
    if m.grid_torque:
      make_grid += "-torque"    
      shm_grid += "-torque"  
    if m.grid_omp:
      make_grid += "-omp"
    for pnr in range(npart):    
      part = "partner-%d" % (pnr+1)
      s = "%s %sr.pdb %s-interior.vol" % (make_grid, part, part)
      s += " $ATTRACTDIR/parmu.par"
      s += " %f %f" % (m.grid_plateau_distance, m.grid_neighbour_distance)
      s += " %s.grid" % part
      if m.dielec == "cdie": s += " --cdie"
      ep = int(10000*m.epsilon)/10000.0
      if ep != 15:
        s += " --epsilon %f" % ep
      if not m.grid_potentials: s += " calc-potentials=0"
      ret += s + "\n"
      ret += shm_grid + " %s.grid %s.gridheader\n" % (part,part)
    ret += "\n"
    
    dock = "python $ATTRACTDIR/protocols/attract.py start.dat"
    dock += " $ATTRACTDIR/parmu.par"
    dock += " " + pdbfiles
    if m.use_grids:
      for n in range(ngrid):
        dock += " --grid %d partner-%d.gridheader" % (n+1, n+1)
        #TODO: homodimers etc.

    if m.gravity:
      dock += " --gravity 1"    #TODO 
    dock += " --rstk %f" % m.rstk
    if m.dielec == "cdie": dock += " --cdie"
    ep = int(10000*m.epsilon)/10000.0
    if ep != 15:
      dock += " --epsilon %f" % ep
    dock += " --jobsize %d" % m.jobsize  
    dock += " --np %d" % m.threads
    dock += " --output dock.dat"
    ret += "\n" + dock + "\n\n"
    
    ret += "$ATTRACTDIR/bin/shm-clean\n\n"
    
  return ret  
"""
