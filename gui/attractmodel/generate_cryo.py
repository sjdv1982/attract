from __future__ import print_function
from math import *

cmd_dock = """python %(protocol)s $inp $dock0 $stage --chunks $threads  --np $threads --output $outp"""

cmd_dockmin = """python %(protocol)s $inp $dock0 $stage --chunks $threads  --np $threads --output $outp-unfiltered
python $ATTRACTDIR/../protocols/attract.py $outp-unfiltered $dock0 $stagemin --chunks $threads  --np $threads --output $outp-min
python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp
"""

script_score = """score() {
  outp=$1
  sorted=$2
  topstruc=$3
  if [ ! -s $topstruc ]; then
%s
      $ATTRACTTOOLS/top $sorted $ntop > $topstruc
  fi
}

"""

script_gvm_noaxsym = """      python $ATTRACTDIR/gvm.py $mapfile $score_threshold $outp $complex | awk '{print "Energy:", $1}' > $outp-gvm1
      python $ATTRACTTOOLS/fill-energies.py $outp $outp-gvm1 > $outp-gvm2
      python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted"""

script_gvm_axsym = """      ./axsym.sh $outp > $outp-axsym
      python $ATTRACTDIR/gvm.py $mapfile $score_threshold $outp-axsym $complex | awk '{print "Energy:", $1}' > $outp-gvm1
      python $ATTRACTTOOLS/fill-energies.py $outp $outp-gvm1 > $outp-gvm2
      python $ATTRACTTOOLS/sort.py $outp-gvm2 --rev > $sorted"""


script_prepare="""prepare(){
  topstruc=$1
  inp=$2
  if [ ! -s $inp ]; then
    python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
    python $ATTRACTDIR/gvm.py $mapfile $score_threshold $topstruc-combined $complex | awk '{print "Energy:", $1}' > $topstruc-gvm1
    python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
    python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc --rev > $topstruc-topcombined
    python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $inp
  fi
}

"""


script_prepare_one_molecule ="""prepare(){
  topstruc=$1
  inp=$2
  if [ ! -s $inp ]; then
    python $ATTRACTTOOLS/monte.py $topstruc $clone > $inp
  fi
}

"""


script_tabufilter="""tabufilter() {
  inp=$1
  ntabu=$2
  outp=$3
  for n in `seq $ntabu`; do
    if [ ! -s $inp-rmsd-tabu$n ]; then
      $ATTRACTDIR/rmsd $inp %s | grep RMSD > $inp-rmsd-tabu$n &
    fi
  done
  wait
  tabus=`seq $ntabu | awk -v inp=$inp '{print inp "-rmsd-tabu" $1}'`
  python $ATTRACTTOOLS/filter-tabu.py $inp %.6f $tabus > $outp-filtered
  $ATTRACTTOOLS/top $outp-filtered $ntop > $outp
}

"""

script_collect =  """collect() {
  inp=$1
  pattern=$2
  $ATTRACTTOOLS/top $inp 1 > best-$pattern.dat
  $ATTRACTTOOLS/top $inp 10 > top10-$pattern.dat
  $ATTRACTDIR/collect best-$pattern.dat %s > best-$pattern.pdb
%s  $ATTRACTDIR/collect top10-$pattern.dat %s > top10-$pattern.pdb
}

"""

script_refine = """refine() {
  run=$1
  inp=$2
  pattern=$3
  ntabu=$4
  for it in `seq $iter`; do
    echo '******************************************************************************************************************************'
    echo Run $run Refinement iteration $it
    echo '******************************************************************************************************************************'

    newinp=$pattern-preit$it.dat
    if [ ! -s $newinp ]; then
      python $ATTRACTTOOLS/monte.py $inp $clone > $newinp
    fi
    outp=$pattern-it$it.dat
    sorted=$pattern-it$it-sorted.dat
    topstruc=$pattern-it$it-topstruc.dat
    if [ ! -s $outp ]; then
      python $ATTRACTDIR/../protocols/attract.py $newinp $dock0 $dock_stage5 --chunks $threads  --np $threads --output $outp
    fi
    if [ ! -s $topstruc ]; then
%s
      if [ $ntabu -eq 0 ]; then
        $ATTRACTTOOLS/top $sorted $ntop > $topstruc
      else
        tabufilter $sorted $ntabu $topstruc
      fi
    fi
    inp=$topstruc
  done
}

"""

script_it_initial = """for it in `seq %d %d`; do
  echo '******************************************************************************************************************************'
  echo Run $run Initial sampling iteration $it '(%s stage, 250 steps)'
  echo '******************************************************************************************************************************'
  outp=dock-run$run-it$it.dat
  sorted=dock-run$run-it$it-sorted.dat
  topstruc=dock-run$run-it$it-topstruc.dat
  nextit=$((it+1))
  if [ ! -s $outp ]; then
   stage=$dock_stage%d
   %s
  fi
  score $outp $sorted $topstruc
  newinp=dock-run$run-preit$nextit.dat
  prepare $topstruc $newinp
  inp=$newinp
done
"""
#(1, 1, first, 1, cmd_dock)
#(2, 4, second, 2, cmd_dock)

script_runs_tabu="""for run in `seq %d %d`; do
  echo '******************************************************************************************************************************'
  echo Run $run
  echo '******************************************************************************************************************************'
  prevrun=$((run-1))
  newinp=dock-run$run-tabufilter.dat
  if [ ! -s $newinp ]; then
    tabufilter %s $prevrun $newinp
  fi
  outp=dock-refine$run-it$iter-topstruc.dat
  if [ ! -s $outp ]; then
    refine $run $newinp dock-refine$run $prevrun
    collect $outp refine$run
  fi
done
"""
#(2, 6, dock-refine1-it1-sorted.dat)
#(7, 10, monocombine-start.dat)

script_recombination = """echo '******************************************************************************************************************************'
echo 'Recombination for run %d-%d'
echo '******************************************************************************************************************************'

if [ ! -s monocombine.dat ]; then
  cat /dev/null > monocombine.list
  for run in `seq %d`; do
    for lig in `seq %d`; do
      f=monocombine-refine$run-lig$lig.dat
      python $ATTRACTTOOLS/monocombine.py $lig dock-run1-it1-sorted.dat best-refine$run.dat > $f
      echo $f >> monocombine.list
    done
  done
  $ATTRACTTOOLS/add `cat monocombine.list` > monocombine.dat
fi
if [ ! -s monocombine-tabu.dat ]; then
  tabufilter monocombine.dat %d monocombine-tabu.dat
  $ATTRACTTOOLS/top monocombine-tabu.dat-filtered $nstruc > monocombine-tabu.dat
fi

if [ ! -s monocombine-start.dat ]; then
      outp=monocombine-tabu-dock.dat
      sorted=monocombine-tabu-dock-sorted.dat
      python $ATTRACTDIR/../protocols/attract.py monocombine-tabu.dat $dock0 $dock_stage5 --chunks $threads  --np $threads --output $outp
%s
      $ATTRACTTOOLS/top monocombine-tabu-dock-sorted.dat $ntop2 > monocombine-start.dat
fi
"""

#(7, 10, 6, 2, 6, script_gvm_noaxsym)
script_main="""#### Main protocol

inp=randsearch.dat

run=1
%(stage1)s
%(stage2)s
%(stage3)s
%(stage4)s
inp=dock-refine1-it$iter-topstruc.dat
if [ ! -s $inp ]; then
  refine 1 dock-run1-it%(initial_stages)d-topstruc.dat dock-refine1 0
fi
collect $inp refine1
%(runs_tabu1)s
%(recombination)s
%(runs_tabu2)s
python $ATTRACTTOOLS/select-em.py %(nruns)d $iter > result-all.dat
$ATTRACTTOOLS/top result-all.dat 100 > result-top100.dat
"""

def generate_cryo_init(m):
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
    ret += "mapmask_threshold=%.6f\n" % m.mapmask_threshold
    ret += "nstruc=%d\n" % m.nstruc
    ret += "ntop=%d\n" % m.ntop
    ret += "ntop2=%d #ntop x 10 (10 clones)\n" % (10 * m.ntop)
    ret += "iter=%d #refinement iterations\n" % m.iterations[4]
    ret += "score_threshold=%.6f\n\n" % m.score_threshold

    filenames = []
    pdbnames = []
    collectnames = []
    tabu_pdbnames = []
    collectnames_final = []
    pdbnames3 = set(("partners.pdb","partners-axsym.pdb"))
    for pnr, p in enumerate(m.partners):
        chain = chr(ord('A') + pnr)
        assert p.code is None #TODO: not implemented
        #TODO: select chain
        pdbname = p.pdbfile.name
        pdbname2 = os.path.split(pdbname)[1]
        pdbname3 = os.path.splitext(pdbname2)[0]
        pdbname3_0 = pdbname3
        pcount = 0
        while pdbname3 in pdbnames3:
            pcount += 1
            pdbname3 = pdbname3_0 + "-" + str(pcount)
        pdbnames3.add(pdbname3)
        pdbname4 = pdbname3 + ".pdb"
        pdbname_heavy = pdbname3 + "-heavy.pdb"
        collectnames.append(pdbname_heavy)
        pdbname_reduced = pdbname3 + "r.pdb"
        if pdbname not in pdbnames:
            if pdbname4 != pdbname:
                ret += "cat %s > %s\n" % (pdbname, pdbname4)
            ret += "python $ATTRACTDIR/../allatom/aareduce.py --heavy --chain %s %s %s > /dev/null\n" % (chain, pdbname4, pdbname_heavy)
            ret += "python $ATTRACTTOOLS/reduce.py --chain %s %s %s > /dev/null\n" % (chain, pdbname_heavy, pdbname_reduced)
        pdbnames.append(pdbname)
        filenames.append(pdbname_reduced)
        tabu_pdbnames.append(pdbname_reduced)
        tabu_pdbnames.append("tabu-refine$n-%s.pdb" % chain)
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
            axstr = " %d %d %.6f %.6f %.6f 0 0 0" % (ax.molecule, ax.fold, ax.axis.x, ax.axis.y, ax.axis.z)
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
        for n in range(len(m.partners)):
            f = filenames[n]
            f2 = collectnames[n]
            ret += "cat %s >> partners-axsym.pdb\n" % f
            ret += "echo TER >> partners-axsym.pdb\n"
            collectnames_final.append(f2)
        for n in range(len(m.axsymmetry)):
            ax = m.axsymmetry[n]
            f = filenames[ax.molecule-1]
            f2 = collectnames[ax.molecule-1]
            copies = ax.fold
            copies_done = 1
            if ax.molecule in axcopies:
                copies *= axcopies[ax.molecule]
                copies_done *= axcopies[ax.molecule]
            axcopies[ax.molecule] = copies
            for nn in range(copies-copies_done):
                collectnames_final.append(f2)
                ret += "cat %s >> partners-axsym.pdb\n" % f
                ret += "echo TER >> partners-axsym.pdb\n"
        ret += "complex=partners-axsym.pdb\n\n"
    else:
        collectnames_final = collectnames
        ret += "complex=partners.pdb\n"

    if m.mapmass is None:
        ret += "mapmass=`python $ATTRACTTOOLS/mass.py $complex`\n"
    else:
        ret += "mapmass=%.6f\n" % m.mapmass
    ret += "python $ATTRACTTOOLS/em/mapsumset-smart.py $mapfile0 $mapfile $mapmass\n"
    ret += "mask1=map-scale-mask1.mask\n"
    ret += "mask2=map-scale-mask2.mask\n"
    ret += "if [ ! -s $mask1 ]; then\n"
    ret += " python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %.6f %.6f $mask1 ${mask1%%%%.*}.sit\n" % (m.mapmask1_voxelsize, m.mapmask1_dimension)
    ret += " python $ATTRACTTOOLS/em/situs2mask.py $mapfile $mapmask_threshold %.6f %.6f $mask2 ${mask2%%%%.*}.sit\n" % (m.mapmask2_voxelsize, m.mapmask2_dimension)
    ret += "fi\n"

    ret += "\n#iteration parameters\n"
    ret += "energy_threshold=%.6f\n" % m.energy_threshold
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

    #stage 1
    dock_stage = "--atomdensitymask \'$mask1\' %.6f --mc --mcscalerot %.6f --mcscalecenter %.6f --mcmax %d" % \
     (m.maskweight[0], m.mcscalerot[0]*fr, m.mcscalecenter[0]*fc, m.mcmax[0])
    if rest:
        dock_stage += rest + " --restweight %.6f" % restweights[0]
    if m.gravity:
        dock_stage += " --gravity 1"
        dock_stage += " --rstk %.6f" % m.rstk
    dock_stages.append(dock_stage)

    #stage 2
    dock_stage = "--atomdensitymask \'$mask2\' %.6f --mc --mcscalerot %.6f --mcscalecenter %.6f --mcmax %d" % \
     (m.maskweight[1], m.mcscalerot[1]*fr, m.mcscalecenter[1]*fc, m.mcmax[1])
    if rest:
        dock_stage += rest + " --restweight %.6f" % restweights[0]
    dock_stages.append(dock_stage)

    #stage 3-4
    for n in range(2,4):
        dock_stage = " --atomdensitymask \'$mask2\' %.6f --mc --mcscalerot %.6f --mcscalecenter %.6f --mcmax %d" % \
         (m.maskweight[n], m.mcscalerot[n]*fr, m.mcscalecenter[n]*fc, m.mcmax[n])
        dock_stage_min = " --atomdensitymask \'$mask2\' %.6f" % m.maskweight[n]
        if m.harmonic_restraints_file:
            t = rest + " --restweight %.6f" % restweights[n]
            dock_stage += t
            dock_stage_min += t
        dock_stages.append(dock_stage)
        dock_stages_min.append(dock_stage_min)


    #stage 5
    dock_stage = "--atomdensitymask \'$mask2\' %s" % m.maskweight[4]
    dock_stages.append(dock_stage)

    ret += "dock_stage1=\'%s\'\n" % dock_stages[0]
    ret += "dock_stage2=\'%s\'\n" % dock_stages[1]
    ret += "dock_stage3=\'%s\'\n" % dock_stages[2]
    ret += "dock_stage4=\'%s\'\n" % dock_stages[3]
    ret += "dock_stage5=\'%s\'\n" % dock_stages[4]
    ret += "dock_stage3_min=\'%s\'\n" % dock_stages_min[0]
    ret += "dock_stage4_min=\'%s\'\n" % dock_stages_min[1]

    clonefac = m.nstruc/m.ntop
    noclone = "\'--ori %.6f --trans %.6f\'" % (m.clone_rot*fr, m.clone_center*fc)
    ret += "noclone=%s\n" % noclone
    clone = "\'--ori %.6f --trans %.6f --clone 10 --keepfirst\'" % (m.clone_rot*fr, m.clone_center*fc)
    ret += "clone=%s\n" % clone
    fast = ""
    if len(m.partners) > 1: fast = " --fast"
    radius = ""
    if m.randsearch_radius != 35: radius = " --radius %.6f" % m.randsearch_radius
    ret += "\nif [ ! -s randsearch.dat ]; then\n python $ATTRACTTOOLS/randsearch.py %d $nstruc%s%s > randsearch.dat\nfi\n"  % (len(m.partners), fast, radius)
    ret += "\n"
    return ret, filenames, tabu_pdbnames, collectnames_final

def generate_cryo(m):
    assert m.score_method == "gvm" #TODO: not implemented
    ret, filenames, tabu_pdbnames, collectnames_final = generate_cryo_init(m)
    if len(m.axsymmetry):
        ret += script_score % script_gvm_axsym
    else:
        ret += script_score % script_gvm_noaxsym
    if len(m.partners) == 1:
      ret += script_prepare_one_molecule
    else:
      ret += script_prepare
    ret += script_tabufilter % (" ".join(tabu_pdbnames), m.tabu_dist)
    collect_str = ""
    for n in range(len(m.partners)):
        ch = chr(ord('A')+n)
        collect_str += "  awk '$1 == \"ATOM\" && substr($0,22,1) == \"%s\"' best-$pattern.pdb > tabu-$pattern-%s.pdb\n" % (ch, ch)
    filenames_str = " ".join(filenames)
    filenames_collect = filenames_str
    if len(m.partners) == 1:
        filenames_collect = "partners.pdb"
    ret += script_collect % (filenames_collect, collect_str, filenames_collect)
    if len(m.axsymmetry):
        script_gvm = script_gvm_axsym
    else:
        script_gvm = script_gvm_noaxsym
    ret += script_refine % script_gvm
    ordinals = "first", "second", "third", "fourth", "fifth"
    cumit = [sum(m.iterations[:n]) for n in range(6)]
    protocol="$ATTRACTDIR/../protocols/attract.py"
    if m.use_gpu:
        protocol="$ATTRACTDIR/../attract-em-gpu/neo-attract-em.py"
    main_params = {}
    main_params["initial_stages"] = cumit[4]
    for n in range(5):
        s = ""
        if m.iterations[n]:
            start = cumit[n] + 1
            end = cumit[n+1]
            cmdd0 = cmd_dockmin if n >=2 else cmd_dock
            cmdd = cmdd0 % {"protocol":protocol}
            s = script_it_initial % (start, end, ordinals[n], n+1, cmdd)
        main_params["stage%d" %(n+1)] = s
    s = ""
    if m.tabu1:
        s = script_runs_tabu % (2, 1 + m.tabu1, "dock-refine1-it1-sorted.dat")
    main_params["runs_tabu1"] = s
    s = ""
    if m.tabu2:
        p = (2 + m.tabu1, 1 + m.tabu1 + m.tabu2, 1 + m.tabu1, len(m.partners),
             1 + m.tabu1, script_gvm)
        s = script_recombination %  p
    main_params["recombination"] = s
    nruns = 1 + m.tabu1 + m.tabu2
    main_params["nruns"] = nruns
    s = ""
    if m.tabu2:
        s = script_runs_tabu % (2 + m.tabu1, nruns, "monocombine-start.dat")
    main_params["runs_tabu2"] = s
    ret += script_main % main_params
    collect_final_str = " ".join(collectnames_final)
    if len(m.axsymmetry):
        ret += "./axsym.sh result-top100.dat > result-top100.dat-axsym\n"
        ret += "$ATTRACTDIR/collect result-top100.dat-axsym %s > result-top100.pdb\n" % collect_final_str
    else:
        ret += "$ATTRACTDIR/collect result-top100.dat %s > result-top100.pdb\n" % collect_final_str
    return ret
