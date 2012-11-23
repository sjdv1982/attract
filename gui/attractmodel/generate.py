#TODO: shm (grids and np > 1)
import os

def generate(m):
  ret = ""
  if (m.header is not None):
    ret = m.header + "\n\n"
  ret += "if [ 1 -eq 1 ]; then ### move and change to disable parts of the protocol\n"

  modes_any = any((p.modes_file for p in m.partners))
  if modes_any:
    ret += """
echo '**************************************************************'
echo 'Assemble modes file...'
echo '**************************************************************'
cat /dev/null > hm-all.dat
"""      
    for p in m.partners:
      if p.nr_modes is None:
        ret += "echo 0 >> hm-all.dat\n"
      elif p.nr_modes == 10:
        ret += "cat %s >> hm-all.dat\n" % p.modes_file.name
      else:
        ret += "#select first %d modes...\n" % (p.nr_modes)
        ret += "awk 'NF == 2 && $1 == %d{exit}{print $0}' %s >> hm-all.dat\n" % \
         (p.nr_modes+1, p.modes_file.name)
    ret += "\n"    
  filenames = []
  reduce_any = False
  for pnr,p in enumerate(m.partners):
    if p.pdb.mode == "download":
      raise Exception("Not implemented")
    if p.is_reduced == False:
      if reduce_any == False:
        ret += """
echo '**************************************************************'
echo 'Reduce partner PDBs...'
echo '**************************************************************'
"""      
        reduce_any = True
      pdbname = p.pdb.pdbfile.name  
      ret += "$ATTRACTDIR/reduce %s >& /dev/null\n" % pdbname
      pdbname_reduced = pdbname[:-4] + "r.pdb"
      filenames.append(pdbname_reduced)
    else:  
      filenames.append(p.pdb.pdbfile.name) 
    #TODO: normal modes  
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
  if m.fix_receptor: params += " --fix-receptor" 
  if modes_any: params += " --modes hm-all.dat"
  params += "\""
  ret += """
#docking parameters
params=%s
"""  % params
  if m.np > 1:
    ret += """
#parallelization parameters
parals="--jobsize %d --np %d"
"""  % (m.jobsize, m.np)
  
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
      ret += "cat %s > rotation.dat\n" % rotfil
      if m.translations_file is not None:
        ret += "cat %s > translate.dat\n" % m.translations_file.name
      else:
        ret += "$ATTRACTDIR/translate %s %s > translate.dat\n" % \
         (filenames[0], filenames[1])
      ret += "$ATTRACTDIR/systsearch > systsearch.dat\n"
      ret += "start=systsearch.dat\n"
    else:
      ret += """
#starting structures
start=%s 
""" % m.start_structures_file.name
  elif m.search == "random":
    ret += "python $ATTRACTTOOLS/randsearch.py %d %d > randsearch.dat\n" % \
     (len(m.partners), m.structures)
    ret += "start=randsearch.dat\n"    
  else:
    raise Exception("Unknown value")
  ret += "\n"  
  ret += """
echo '**************************************************************'
echo 'Docking'
echo '**************************************************************'
"""
  inp = "$start"
  ordinals = ["1st", "2nd", "3rd",] + ["%dth" % n for n in range(4,51)]
  iterations = []
  if m.iterations is None or len(m.iterations) == 0:
    iterations.append([None, None, False, False])
  else:
    for it in m.iterations:
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
    if rcut is not None: itparams += " --rcut %s" % str(rcut)
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
    ret += "%s %s $params %s %s %s\n" % (attract, inp, itparams, tail, outp)
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
    ignorens = ""
    any_ensemble = any((p.ensemble_size is not None for p in m.partners))
    if any_ensemble and m.deredundant_ignorens: ignorens = " --ignorens"
    ret += "python $ATTRACTTOOLS/deredundant.py %s%s > %s\n" % (result, ignorens, outp)
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
          else:
            raise Exception("TODO: ensembles")  
          result = outp
          outp, outp2 = outp2, outp
  if m.calc_lrmsd:
    #TODO: ensembles
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
        filename = p.rmsd_pdb
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
    #TODO: modes, ens
    lrmsd_allfilenames = []
    for f1, f2 in zip(lrmsd_filenames, lrmsd_refenames):
      lrmsd_allfilenames.append(f1)
      lrmsd_allfilenames.append(f2)
    lrmsd_allfilenames = " ".join(lrmsd_allfilenames)
    lrmsdresult = os.path.splitext(result0)[0] + ".lrmsd"
    par = ""
    if modes_any and not deflex_any: par = " --modes hm-all.dat"
    ret += "$ATTRACTDIR/lrmsd %s %s%s > %s\n" % (result, lrmsd_allfilenames, par, lrmsdresult)
    ret += "\n"

  if m.calc_lrmsd or m.calc_irmsd:
    if deflex_any:
      ret += "rm -f $tmpf $tmpf2\n"

  if m.collect:
    #TODO: ensembles, modes
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
    if modes_any: par = " --modes hm-all.dat"    
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
