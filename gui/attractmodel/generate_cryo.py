from __future__ import print_function
from math import *

itscript_axsym = """
outp=dock-it%(it)s.dat
sorted=dock-it%(it)s-sorted.dat
topstruc=dock-it%(it)s-topstruc.dat
python $ATTRACTDIR/../protocols/attract.py $inp $dock0 %(dock)s --chunks $threads1  --np $threads1 --output $outp
python $ATTRACTTOOLS/filter-energy.py $outp $energy_threshold > $outp-filtered
./axsym.sh $outp-filtered > $outp-axsym
python $ATTRACTDIR/gvm.py $mapfile $score_threshold $outp-axsym $complex | awk '{print "Energy:", $1}' > $outp-gvm1
python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted
$ATTRACTTOOLS/top $sorted $ntop > $topstruc
""" 

itscript_noaxsym = """
outp=dock-it%(it)s.dat
sorted=dock-it%(it)s-sorted.dat
topstruc=dock-it%(it)s-topstruc.dat
python $ATTRACTDIR/../protocols/attract.py $inp $dock0 %(dock)s --chunks $threads1  --np $threads1 --output $outp
python $ATTRACTTOOLS/filter-energy.py $outp $energy_threshold > $outp-filtered
python $ATTRACTDIR/gvm.py $mapfile $score_threshold $outp-filtered $complex | awk '{print "Energy:", $1}' > $outp-gvm1
python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted
$ATTRACTTOOLS/top $sorted $ntop > $topstruc
""" 

def generate_cryo(m):
  import os
  ret = """#!/bin/bash
set -u -e
"""
  np = len(m.partners)
  ret += "threads1=%d\n" % m.threads1
  ret += "threads2=%d\n\n" % m.threads2
  ret += "name=%s\n" % m.runname
  ret += "nbodies=%d\n" % np
  mapfilename = m.mapfilename
  if mapfilename is None: 
    mapfilename = m.mapfile.name
  ret += "mapfile=%s\n" % mapfilename
  ret += "mapmask_threshold=%s\n" % m.mapmask_threshold
  ret += "nstruc=%d\n" % m.nstruc
  ret += "ntop=%d\n" % m.ntop
  ret += "score_threshold=%s\n\n" % m.score_threshold
  
  filenames = []
  pdbnames = []
  pdbnames3 = set(("partners.pdb","partners-axsym.pdb"))
  for pnr,p in enumerate(m.partners):
    assert p.code is None #TODO: not implemented
    #TODO: select chain
    pdbname = p.pdbfile.name
    pdbname2 = os.path.split(pdbname)[1]
    if pdbname not in pdbnames: 
      if pdbname2.count(".") > 1:
        raise ValueError("The 'reduce' program does not support PDB files with double extension: '%s'" % pdbname2)
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
      pdbname_reduced = pdbname3 + "r.pdb"
      ret += "$ATTRACTDIR/reduce %s > /dev/null\n" % pdbname4            
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
  
  if m.mapmass is None:
    ret += "mapmass=`python $ATTRACTTOOLS/mass.py $complex`\n"
  else:
    ret += "mapmass=%s" % m.mapmass
  ret += "mask1=map-scale-mask1.mask\n"
  ret += "mask2=map-scale-mask2.mask\n"
  ret += "python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %s %s $mask1 ${mask1%%%%.*}.sit\n" % (m.mapmask1_voxelsize, m.mapmask1_dimension)
  ret += "python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %s %s $mask2 ${mask2%%%%.*}.sit\n" % (m.mapmask2_voxelsize, m.mapmask2_dimension)
  
  ret += "\n#iteration parameters\n"
  ret += "energy_threshold=%s\n" % m.energy_threshold
  dock0 = "$ATTRACTDIR'/../parmw.par partners.pdb --ghost" + axpar + "'"
  assert m.restraints_file is None #TODO
  ret += "dock0=%s\n" % dock0 
  dock_stages = []
  fr = m.global_scale_rot
  fc = m.global_scale_trans
  dock_stage = "--atomdensitymask \'$mask1\' %s --mc --mcscalerot %s --mcscalecenter %s --mcmax %s" % \
   (m.maskweight[0], m.mcscalerot[0]*fr, m.mcscalecenter[0]*fc, m.mcmax[0])
  if m.gravity:
    dock_stage += " --gravity 1"   
  if m.gravity or m.restraints_file:
    dock_stage += " --rstk %s" % m.rstk    
  dock_stages.append(dock_stage) 
  for n in range(1,4):
    dock_stage = " --atomdensitymask \'$mask2\' %s --mc --mcscalerot %s --mcscalecenter %s --mcmax %s" % \
     (m.maskweight[n], m.mcscalerot[n]*fr, m.mcscalecenter[n]*fc, m.mcmax[n])
    dock_stages.append(dock_stage)
  dock_stage = "--atomdensitymask \'$mask2\' %s" % m.maskweight[4]
  dock_stages.append(dock_stage)
  ret += "dock_stage1=\'%s\'\n" % dock_stages[0]
  ret += "dock_stage2=\'%s\'\n" % dock_stages[1]
  ret += "dock_stage3=\'%s\'\n" % dock_stages[2]
  ret += "dock_stage4=\'%s\'\n" % dock_stages[3]
  ret += "dock_stage5=\'%s\'\n" % dock_stages[4]
  clonefac = m.nstruc/m.ntop
  clone1 = "\'--ori %s --trans %s --clone %d --keepfirst\'" % (m.clone_rot[0]*fr, m.clone_center[0]*fc, clonefac)
  clone2 = "\'--ori %s --trans %s --clone %d --keepfirst\'" % (m.clone_rot[1]*fr, m.clone_center[1]*fc, clonefac)
  ret += "clone1=%s\n" % clone1
  ret += "clone2=%s\n" % clone2
  fast = ""
  if len(m.partners) > 1: fast = " --fast"
  ret += "\npython $ATTRACTTOOLS/randsearch.py %d %d%s > randsearch.dat\n"  % (len(m.partners), m.nstruc, fast)
  ret += "inp=randsearch.dat\n\n"

  ret += "#iterations\n"
  if len(m.axsymmetry):
    itscript = itscript_axsym
  else:
    itscript = itscript_noaxsym
  for it in range(1,m.iterations+1):
    dock = "$dock_stage1"
    clone = "$clone1"
    if it > 1:      
      if it <= m.mc_iterations:
        if it == 2:
          dock = "$dock_stage2"
        else:  
          stage34 = m.mc_iterations - 2
          if (it-3) < stage34/2.0:
            dock = "$dock_stage3"
          else:
            dock = "$dock_stage4"
        if it == m.mc_iterations: clone = "$clone2"    
      else:
        dock = "$dock_stage5"
        clone = "$clone2"
    ret += itscript % {"dock": dock, "it": it}    
    assert m.score_method == "gvm" #TODO: non-GVM
    topstrucf = "$topstruc"
    if len(m.partners) > 1 and it <= m.mc_iterations:
      ret += "python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined\n"
      ret += "python $ATTRACTDIR/gvm.py $mapfile $score_threshold $topstruc-combined $complex | awk '{print \"Energy:\", $1}' > $topstruc-gvm1\n"
      ret += "python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2\n"
      ret += "python $ATTRACTTOOLS/sort.py $topstruc-gvm2 --rev > $topstruc-sorted\n"
      ret += "$ATTRACTTOOLS/top $topstruc-sorted $ntop > $topstruc-recombined\n"      
      topstrucf = "$topstruc-recombined"
    if it < m.iterations:
      ret += "newinp=dock-preit%s.dat\n" % (it+1)
      ret += "python $ATTRACTTOOLS/monte.py %s %s > $newinp\n" % (topstrucf, clone)
      ret += "inp=$newinp\n"
  return ret