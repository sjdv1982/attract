from __future__ import print_function
from math import *


itscript1 = """outp=dock-run%(run)s-it%(it)s.dat
sorted=dock-run%(run)s-it%(it)s-sorted.dat
topstruc=dock-run%(run)s-it%(it)s-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python $ATTRACTDIR/../protocols/attract.py $inp $dock0 %(dock)s --chunks $threads  --np $threads --output $outp
 python $ATTRACTTOOLS/filter-energy.py $outp $energy_threshold > $outp-filtered
"""

itscript1min = """outp=dock-run%(run)s-it%(it)s.dat
sorted=dock-run%(run)s-it%(it)s-sorted.dat
topstruc=dock-run%(run)s-it%(it)s-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python $ATTRACTDIR/../protocols/attract.py $inp $dock0 %(dock)s --chunks $threads  --np $threads --output $outp
 python $ATTRACTDIR/../protocols/attract.py $outp $dock0 %(dock_min)s --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
"""

itscript1_dummy = """outp=dock-run%(run)s-it%(it)s.dat
sorted=dock-run%(run)s-it%(it)s-sorted.dat
topstruc=dock-run%(run)s-it%(it)s-topstruc.dat
if [ ! -s $outp-filtered ]; then
 ln -s dock-run1-it%(it)s.dat-filtered $outp-filtered 
"""

itscript1_dummy2 = """outp=dock-run%(run)s-it%(it)s.dat
sorted=dock-run%(run)s-it%(it)s-sorted.dat
topstruc=dock-run%(run)s-it%(it)s-topstruc.dat
if [ ! -s $outp-filtered ]; then
 ln -s dock-run1-it%(it)s-topstruc.dat-topcombined $outp-filtered 
"""

itscript2_axsym = """if [ ! -s $topstruc ]; then
 ./axsym.sh %s > $outp-axsym
 python $ATTRACTDIR/gvm.py $mapfile $score_threshold $outp-axsym $complex | awk '{print "Energy:", $1}' > $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py %s $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted
 $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
""" 

itscript2_noaxsym = """if [ ! -s $topstruc ]; then
 python $ATTRACTDIR/gvm.py $mapfile $score_threshold %s $complex | awk '{print "Energy:", $1}' > $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py %s $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted
 $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
""" 

def generate_cryo(m):
  import os
  ret = """#!/bin/bash -i
set -u -e
"""
  np = len(m.partners)
  ret += "threads=%d\n" % m.threads
  ret += "name=%s\n" % m.runname
  ret += "nbodies=%d\n" % np
  mapfilename = m.mapfilename
  if mapfilename is None: 
    mapfilename = m.mapfile.name
  mapfilename2 = "map-scaled.sit"
  if mapfilename == mapfilename2:
    mapfilename2 = "map-scaled-scaled.sit"
  ret += "mapfile0=%s\n" % mapfilename
  ret += "mapfile=%s\n" % mapfilename2
  ret += "mapmask_threshold=%s\n" % m.mapmask_threshold
  ret += "nstruc=%d\n" % m.nstruc
  ret += "ntop=%d\n" % m.ntop
  ret += "score_threshold=%s\n\n" % m.score_threshold
  
  filenames = []
  pdbnames = []
  collectnames = []
  pdbnames3 = set(("partners.pdb","partners-axsym.pdb"))
  for pnr, p in enumerate(m.partners):
    chain = chr(ord('A') + pnr)
    assert p.code is None #TODO: not implemented
    #TODO: select chain
    pdbname = p.pdbfile.name
    pdbname2 = os.path.split(pdbname)[1]
    if pdbname not in pdbnames: 
      pdbname3 = os.path.splitext(pdbname2)[0]
      pdbname3_0 = pdbname3
      pcount = 0
      while pdbname3 in pdbnames3:
        pcount += 1
        pdbname3 = pdbname3_0 + "-" + str(pcount)
      pdbnames3.add(pdbname3)
      pdbname4 = pdbname3 + ".pdb"
      if pdbname4 != pdbname:
        ret += "cat %s > %s\n" % (pdbname, pdbname4)          
      collectnames.append(pdbname4)
      pdbname_heavy = pdbname3 + "-heavy.pdb"
      pdbname_reduced = pdbname3 + "r.pdb"
      ret += "python $ATTRACTDIR/../allatom/aareduce.py --heavy --chain %s %s %s > /dev/null\n" % (chain, pdbname4, pdbname_heavy)
      ret += "python $ATTRACTTOOLS/reduce.py --chain %s %s %s > /dev/null\n" % (chain, pdbname_heavy, pdbname_reduced)
      pdbnames.append(pdbname)
    else:
      pdbname_reduced = filenames[pdbnames.index(pdbname)]
    filenames.append(pdbname_reduced)
  ret += "echo %d > partners.pdb\n" % np
  for f in filenames:
    ret += "cat %s >> partners.pdb\n" % f
    ret += "echo TER >> partners.pdb\n"
  ret += "\n"
  
  structure_partners = open("partners.dat", "w")
  print("#pivot auto", file=structure_partners)
  print("#1", file=structure_partners)
  for n in range(np):
    print("0 0 0 0 0 0", file=structure_partners)
  structure_partners.close()
  axpar = ""
  if len(m.axsymmetry):
    axsym = open("axsym.sh", "w")
    axsymstr = "$ATTRACTDIR/axsymmetry $1 %d" % np
    npax = np
    axcopies = {}
    for n in range(len(m.axsymmetry)):
      ax = m.axsymmetry[n]
      axstr = " %d %d %s %s %s 0 0 0" % (ax.molecule, ax.fold, ax.axis.x, ax.axis.y, ax.axis.z)
      axsymstr += axstr
      axpar += " --axsym" + axstr
      copies = ax.fold
      copies_done = 1
      if ax.molecule in axcopies: 
        copies *= axcopies[ax.molecule]
        copies_done *= axcopies[ax.molecule]
      axcopies[ax.molecule] = copies        
      npax += copies - copies_done
    print(axsymstr, file=axsym)
    axsym.close()    
    os.system("chmod +x axsym.sh")    
    
    axcopies = {}
    ret += "echo %d > partners-axsym.pdb\n" % npax
    for f in filenames:
      ret += "cat %s >> partners-axsym.pdb\n" % f
      ret += "echo TER >> partners-axsym.pdb\n"    
    for n in range(len(m.axsymmetry)):
      ax = m.axsymmetry[n]
      f = filenames[ax.molecule-1]
      copies = ax.fold
      copies_done = 1
      if ax.molecule in axcopies: 
        copies *= axcopies[ax.molecule]
        copies_done *= axcopies[ax.molecule]
      axcopies[ax.molecule] = copies  
      for nn in range(copies-copies_done):
        ret += "cat %s >> partners-axsym.pdb\n" % f
        ret += "echo TER >> partners-axsym.pdb\n"    
    ret += "complex=partners-axsym.pdb\n\n"
  else:
    ret += "complex=partners.pdb\n"     

  tabunames = []
  if m.mapmass is None:
    ret += "mapmass=`python $ATTRACTTOOLS/mass.py $complex`\n"
  else:
    ret += "mapmass=%s\n" % m.mapmass
  ret += "python $ATTRACTTOOLS/em/mapsumset-smart.py $mapfile0 $mapfile $mapmass\n" 
  ret += "mask1=map-scale-mask1.mask\n"
  ret += "mask2=map-scale-mask2.mask\n"
  ret += "if [ ! -s $mask1 ]; then\n"
  ret += " python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %s %s $mask1 ${mask1%%%%.*}.sit\n" % (m.mapmask1_voxelsize, m.mapmask1_dimension)
  ret += " python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %s %s $mask2 ${mask2%%%%.*}.sit\n" % (m.mapmask2_voxelsize, m.mapmask2_dimension)
  ret += "fi\n"
  
  ret += "\n#iteration parameters\n"
  ret += "energy_threshold=%s\n" % m.energy_threshold
  dock0 = "$ATTRACTDIR'/../attract.par partners.pdb --ghost" + axpar + "'"
  rest = ""
  restweights = 0.1, 0.2, 0.6, 0.6, 1.0 
  if m.harmonic_restraints_file is not None:
    raise NotImplementedError #TODO: convert with tbl2attract    
    restfile = "harmonic-restraints.rest"
    rest = " --rest %s" % restfile
  ret += "dock0=%s\n" % dock0 
  dock_stages = []
  dock_stages_min = []
  fr = m.global_scale_rot
  fc = m.global_scale_trans    
  dock_stage = "--atomdensitymask \'$mask1\' %s --mc --mcscalerot %s --mcscalecenter %s --mcmax %s" % \
   (m.maskweight[0], m.mcscalerot[0]*fr, m.mcscalecenter[0]*fc, m.mcmax[0])
  if rest:
    dock_stage += rest + " --restweight %s" % restweights[0]
  if m.gravity:
    dock_stage += " --gravity 1"   
    dock_stage += " --rstk %s" % m.rstk    
  dock_stages.append(dock_stage) 
  for n in range(1,4):
    dock_stage = " --atomdensitymask \'$mask2\' %s --mc --mcscalerot %s --mcscalecenter %s --mcmax %s" % \
     (m.maskweight[n], m.mcscalerot[n]*fr, m.mcscalecenter[n]*fc, m.mcmax[n])
    dock_stage_min = " --atomdensitymask \'$mask2\' %s" % m.maskweight[n]
    if m.harmonic_restraints_file:
      t = rest + " --restweight %s" % restweights[n]
      dock_stage += t
      dock_stage_min += t
    dock_stages.append(dock_stage)
    dock_stages_min.append(dock_stage_min)
  dock_stage = "--atomdensitymask \'$mask2\' %s" % m.maskweight[4]
  dock_stages.append(dock_stage)
  ret += "dock_stage1=\'%s\'\n" % dock_stages[0]
  ret += "dock_stage2=\'%s\'\n" % dock_stages[1]
  ret += "dock_stage3=\'%s\'\n" % dock_stages[2]
  ret += "dock_stage4=\'%s\'\n" % dock_stages[3]
  ret += "dock_stage5=\'%s\'\n" % dock_stages[4]
  if len(m.partners) > 1:
    ret += "dock_stage2_min=\'%s\'\n" % dock_stages_min[0]
    ret += "dock_stage3_min=\'%s\'\n" % dock_stages_min[1]
    ret += "dock_stage4_min=\'%s\'\n" % dock_stages_min[2]
  
  clonefac = m.nstruc/m.ntop
  if len(m.partners) == 1:
    clone1 = "\'--ori %s --trans %s --clone %d --keepfirst\'" % (m.clone_rot[0]*fr, m.clone_center[0]*fc, clonefac)
    ret += "clone1=%s\n" % clone1
  else:
    noclone = "\'--ori %s --trans %s\'" % (m.clone_rot[1]*fr, m.clone_center[1]*fc)
    ret += "noclone=%s\n" % noclone
  clone2 = "\'--ori %s --trans %s --clone %d --keepfirst\'" % (m.clone_rot[1]*fr, m.clone_center[1]*fc, clonefac)  
  ret += "clone2=%s\n" % clone2
  fast = ""
  if len(m.partners) > 1: fast = " --fast"
  radius = ""
  if m.randsearch_radius != 35: radius = " --radius %s" % m.randsearch_radius  
  ret += "\nif [ ! -s randsearch.dat ]; then\n python $ATTRACTTOOLS/randsearch.py %d $nstruc%s%s > randsearch.dat\nfi\n"  % (len(m.partners), fast, radius)
  ret += "inp=randsearch.dat\n\n"

  runs = m.tabu + 1
  if len(m.partners) > 1:
    runs += 6
  for run in range(1, runs+1):
    ret += "#run %d\n" % run  
    if len(m.axsymmetry):
      itscript2 = itscript2_axsym
    else:
      itscript2 = itscript2_noaxsym
    totit = sum(m.iterations)
    cumsumit = []
    cumsumit0 = 0
    for n in range(len(m.iterations)):
      cumsumit0 += m.iterations[n]
      cumsumit.append(cumsumit0)
    for it in range(1,totit+1):
      dock = "$dock_stage1"
      clone = "$clone1"
      dock_min = None
      if run > m.tabu+1:
        dock = "$dock_stage5"
        clone = "$clone2"        
      else:  
        if it > cumsumit[0]:      
          if it <= cumsumit[3]:
            if it <= cumsumit[1]:
              dock = "$dock_stage2"
              dock_min = "$dock_stage2_min"
            elif it <= cumsumit[2]:
              dock = "$dock_stage3"
              dock_min = "$dock_stage3_min"
            else:  
              dock = "$dock_stage4"
              dock_min = "$dock_stage4_min"
            if it == cumsumit[3]: clone = "$clone2"    
          else:
            dock = "$dock_stage5"
            clone = "$clone2"
      if run > 1 and it < cumsumit[1]: 
          continue
      if run > 1 and it == cumsumit[1]:
        if run > m.tabu + 1:
          ret += itscript1_dummy2 % {"run": run, "it" : it}    
        else:  
          ret += itscript1_dummy % {"run": run, "it" : it}    
      elif dock_min is not None:
        ret += itscript1min % {"run": run, "dock": dock, "dock_min" : dock_min, "it": it}            
      else:
        ret += itscript1 % {"run": run, "dock": dock, "it": it}    
      filtered = "$outp-filtered"    
      if len(tabunames):
        rmsdfiles = []
        for n in range(len(tabunames)):
          tfiles = tabunames[n]
          rmsdpdbfiles = []
          for a,b in zip(collectnames, tfiles): rmsdpdbfiles += [a,b]
          rmsdfile = filtered + "-rmsd-tabu%d" % (n+1)
          rmsdfiles.append(rmsdfile)
          ret += " $ATTRACTDIR/rmsd %s %s | grep RMSD > %s\n" % (filtered, " ".join(rmsdpdbfiles), rmsdfile)
        newfiltered = filtered + "-tabu"
        ret += " python $ATTRACTTOOLS/filter-tabu.py %s 10.0 %s > %s\n" % (filtered, " ".join(rmsdfiles), newfiltered)
        filtered = newfiltered      
      ret += "fi\n"  
      ret += itscript2 % (filtered, filtered)
      assert m.score_method == "gvm" #TODO: non-GVM    
      if it < totit:
        ret += "newinp=dock-run%s-preit%s.dat\n" % (run, it+1)
        ret += "if [ ! -s $newinp ]; then\n"
        if len(m.partners) > 1 and run <= (m.tabu+1) and it <= cumsumit[3]:
          ret += " python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined\n"                    
          ret += " python $ATTRACTDIR/gvm.py $mapfile $score_threshold $topstruc-combined $complex | awk '{print \"Energy:\", $1}' > $topstruc-gvm1\n"
          ret += " python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2\n"
          ret += " python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc --rev > $topstruc-topcombined\n"
          ret += " python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp\n"
        else:      
          ret += " python $ATTRACTTOOLS/monte.py $topstruc %s > $newinp\n" %  clone
        ret += "fi\n"
        ret += "inp=$newinp\n\n"
      else:
        ret += "\nif [ ! -s result-run%d.dat ]; then\n ln -s $topstruc result-run%d.dat\nfi\n" % (run, run)
        ret += "$ATTRACTTOOLS/top result-run%d.dat 1 > best-result-run%d.dat\n" % (run, run)
        ret += "$ATTRACTTOOLS/top result-run%d.dat 10 > top10-result-run%d.dat\n" % (run, run)
        collectnamestr = " ".join(collectnames)
        ret += "$ATTRACTDIR/collect best-result-run%d.dat %s > best-result-run%d.pdb\n" % (run, collectnamestr, run)
        t = []
        for n in range(len(m.partners)):
          ch = chr(ord('A')+n)
          f = "tabu-run%d-%s.pdb" % (run, ch)
          ret += "awk '$1 == \"ATOM\" && substr($0,22,1) == \"%s\"' best-result-run%d.pdb > %s\n" % (ch, run, f)
          t.append(f)
        tabunames.append(t)  
        ret += "$ATTRACTDIR/collect top10-result-run%d.dat %s > top10-result-run%d.pdb\n" % (run, collectnamestr, run)
        ret += "\n"
  return ret