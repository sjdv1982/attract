#!/bin/bash
set -u -e
threads=2
name=attract-em
nbodies=2
mapfile=beadmodel.pdb
nstruc=6000
ntop=300
score_threshold=10.0

cp ubA-aa.pdb proteinA.pdb
cp ubB-aa.pdb proteinB.pdb
grep CA -w proteinB.pdb > proteinB-CA.pdb
python $ATTRACTTOOLS/reduce.py proteinA.pdb > /dev/null
python $ATTRACTTOOLS/reduce.py proteinB.pdb > /dev/null
echo 2 > partners.pdb
cat proteinAr.pdb >> partners.pdb
echo TER >> partners.pdb
cat proteinBr.pdb >> partners.pdb
echo TER >> partners.pdb

complex=partners.pdb

mask1=map-scale-mask1.mask
mask2=map-scale-mask2.mask
if [ ! -s $mask1 ]; then
 python $ATTRACTTOOLS/saxs/pdb2mask.py $mapfile 10.0 100.0 $mask1 ${mask1%%.*}.sit
fi
if [ ! -s $mask2 ]; then
 python $ATTRACTTOOLS/saxs/pdb2mask.py $mapfile 5.0 50.0 $mask2 ${mask2%%.*}.sit
fi

#iteration parameters
energy_threshold=10000.0
dock0=$ATTRACTDIR'/../attract.par partners.pdb --ghost'
dock_stage1='--atomdensitymask '$mask1' 8000.0 --mc --mcscalerot 0.19625 --mcscalecenter 5.0 --mcmax 500 --gravity 1 --rstk 0.1'
dock_stage2=' --atomdensitymask '$mask1' 8000.0 --mc --mcscalerot 0.098125 --mcscalecenter 2.5 --mcmax 500'
dock_stage3=' --atomdensitymask '$mask2' 600.0 --mc --mcscalerot 0.03925 --mcscalecenter 1.0 --mcmax 250'
dock_stage4=' --atomdensitymask '$mask2' 1000.0 --mc --mcscalerot 0.03925 --mcscalecenter 1.0 --mcmax 250'
dock_stage5='--atomdensitymask '$mask2' 1000.0'
dock_stage2_min=' --atomdensitymask '$mask1' 8000.0'
dock_stage3_min=' --atomdensitymask '$mask2' 600.0'
dock_stage4_min=' --atomdensitymask '$mask2' 1000.0'
noclone='--ori 0.03925 --trans 1.0'
clone2='--ori 0.03925 --trans 1.0 --clone 20 --keepfirst'
clone3='--ori 0.03925 --trans 1.0 --clone 5 --keepfirst'
if [ ! -s randsearch.dat ]; then
 python $ATTRACTTOOLS/randsearch.py 2 6000 --fast > randsearch.dat
fi
inp=randsearch.dat

#run 1
outp=dock-run1-it1.dat
sorted=dock-run1-it1-sorted.dat
topstruc=dock-run1-it1-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage1 --chunks $threads  --np $threads --output $outp
 python $ATTRACTTOOLS/filter-energy.py $outp $energy_threshold > $outp-filtered
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted
  $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit2.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it2.dat
sorted=dock-run1-it2-sorted.dat
topstruc=dock-run1-it2-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage2 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage2_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2 
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted
  $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit3.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it3.dat
sorted=dock-run1-it3-sorted.dat
topstruc=dock-run1-it3-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage2 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage2_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2 
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted
  $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit4.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it4.dat
sorted=dock-run1-it4-sorted.dat
topstruc=dock-run1-it4-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage3 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage3_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted
  $ATTRACTTOOLS/top $sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit5.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it5.dat
sorted=dock-run1-it5-sorted.dat
topstruc=dock-run1-it5-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage3 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage3_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre  
 #select first cluster only necessary in first iteration
 $ATTRACTTOOLS/top $sorted 1 > taboo.dat
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre
 python $ATTRACTTOOLS/saxs/savechi.py taboo.dat taboochi.npy
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit6.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it6.dat
sorted=dock-run1-it6-sorted.dat
topstruc=dock-run1-it6-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage3 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage3_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
 
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit7.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it7.dat
sorted=dock-run1-it7-sorted.dat
topstruc=dock-run1-it7-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
 
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit8.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it8.dat
sorted=dock-run1-it8-sorted.dat
topstruc=dock-run1-it8-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
 
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit9.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it9.dat
sorted=dock-run1-it9-sorted.dat
topstruc=dock-run1-it9-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
 
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit10.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it10.dat
sorted=dock-run1-it10-sorted.dat
topstruc=dock-run1-it10-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit11.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp

fi
inp=$newinp

outp=dock-run1-it11.dat
sorted=dock-run1-it11-sorted.dat
topstruc=dock-run1-it11-topstruc.dat
if [ ! -s $outp-filtered ]; then
 
 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit12.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it12.dat
sorted=dock-run1-it12-sorted.dat
topstruc=dock-run1-it12-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit13.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it13.dat
sorted=dock-run1-it13-sorted.dat
topstruc=dock-run1-it13-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit14.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it14.dat
sorted=dock-run1-it14-sorted.dat
topstruc=dock-run1-it14-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit15.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it15.dat
sorted=dock-run1-it15-sorted.dat
topstruc=dock-run1-it15-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit16.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it16.dat
sorted=dock-run1-it16-sorted.dat
topstruc=dock-run1-it16-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit17.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it17.dat
sorted=dock-run1-it17-sorted.dat
topstruc=dock-run1-it17-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit18.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it18.dat
sorted=dock-run1-it18-sorted.dat
topstruc=dock-run1-it18-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit19.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it19.dat
sorted=dock-run1-it19-sorted.dat
topstruc=dock-run1-it19-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit20.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it20.dat
sorted=dock-run1-it20-sorted.dat
topstruc=dock-run1-it20-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit21.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it21.dat
sorted=dock-run1-it21-sorted.dat
topstruc=dock-run1-it21-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit22.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it22.dat
sorted=dock-run1-it22-sorted.dat
topstruc=dock-run1-it22-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit23.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it23.dat
sorted=dock-run1-it23-sorted.dat
topstruc=dock-run1-it23-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit24.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it24.dat
sorted=dock-run1-it24-sorted.dat
topstruc=dock-run1-it24-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit25.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it25.dat
sorted=dock-run1-it25-sorted.dat
topstruc=dock-run1-it25-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit26.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it26.dat
sorted=dock-run1-it26-sorted.dat
topstruc=dock-run1-it26-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit27.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it27.dat
sorted=dock-run1-it27-sorted.dat
topstruc=dock-run1-it27-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit28.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it28.dat
sorted=dock-run1-it28-sorted.dat
topstruc=dock-run1-it28-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit29.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it29.dat
sorted=dock-run1-it29-sorted.dat
topstruc=dock-run1-it29-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit30.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it30.dat
sorted=dock-run1-it30-sorted.dat
topstruc=dock-run1-it30-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit31.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it31.dat
sorted=dock-run1-it31-sorted.dat
topstruc=dock-run1-it31-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit32.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it32.dat
sorted=dock-run1-it32-sorted.dat
topstruc=dock-run1-it32-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit33.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it33.dat
sorted=dock-run1-it33-sorted.dat
topstruc=dock-run1-it33-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit34.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it34.dat
sorted=dock-run1-it34-sorted.dat
topstruc=dock-run1-it34-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit35.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it35.dat
sorted=dock-run1-it35-sorted.dat
topstruc=dock-run1-it35-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit36.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it36.dat
sorted=dock-run1-it36-sorted.dat
topstruc=dock-run1-it36-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit37.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp


outp=dock-run1-it37.dat
sorted=dock-run1-it37-sorted.dat
topstruc=dock-run1-it37-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit38.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it38.dat
sorted=dock-run1-it38-sorted.dat
topstruc=dock-run1-it38-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit39.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it39.dat
sorted=dock-run1-it39-sorted.dat
topstruc=dock-run1-it39-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit40.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it40.dat
sorted=dock-run1-it40-sorted.dat
topstruc=dock-run1-it40-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit41.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it41.dat
sorted=dock-run1-it41-sorted.dat
topstruc=dock-run1-it41-topstruc.dat
if [ ! -s $outp-filtered ]; then

 
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi 
newinp=dock-run1-preit42.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it42.dat
sorted=dock-run1-it42-sorted.dat
topstruc=dock-run1-it42-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit43.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it43.dat
sorted=dock-run1-it43-sorted.dat
topstruc=dock-run1-it43-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit44.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it44.dat
sorted=dock-run1-it44-sorted.dat
topstruc=dock-run1-it44-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit45.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it46.dat
sorted=dock-run1-it46-sorted.dat
topstruc=dock-run1-it46-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-1
 python ../neo-attract-em.py $outp-1 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-1
 
 by=25000
 
 
 python $ATTRACTTOOLS/split.py taboo.dat tmp123 8
 
 python ../neo-attract-em.py tmp123-1 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-2 
  if [ -s tmp123-2 ]; then  
 python ../neo-attract-em.py tmp123-2 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-3 
 python ../neo-attract-em.py tmp123-3 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-4 
 python ../neo-attract-em.py tmp123-4 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-5 
 python ../neo-attract-em.py tmp123-5 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-6 
 python ../neo-attract-em.py tmp123-6 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-7 
 python ../neo-attract-em.py tmp123-7 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-8 
 python ../neo-attract-em.py tmp123-8 $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp-9 
 python ../neo-attract-em.py $outp-3 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-3 
 python ../neo-attract-em.py $outp-4 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-4 
 python ../neo-attract-em.py $outp-5 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-5 
 python ../neo-attract-em.py $outp-6 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-6 
 python ../neo-attract-em.py $outp-7 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-7 
 python ../neo-attract-em.py $outp-8 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-8 
  python ../neo-attract-em.py $outp-9 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-9 
 fi
 python ../neo-attract-em.py $outp-2 $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min-2
 
  if [ -s tmp123-2 ]; then  
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 $outp-min-3 $outp-min-4 $outp-min-5 $outp-min-6 $outp-min-7 $outp-min-8 $outp-min-9 > $outp-min 
 fi 
 if  [ ! -s tmp123-2 ]; then 
 $ATTRACTTOOLS/add $outp-min-1 $outp-min-2 > $outp-min 
 fi
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit47.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it47.dat
sorted=dock-run1-it47-sorted.dat
topstruc=dock-run1-it47-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit48.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it48.dat
sorted=dock-run1-it48-sorted.dat
topstruc=dock-run1-it48-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit49.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it49.dat
sorted=dock-run1-it49-sorted.dat
topstruc=dock-run1-it49-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi
newinp=dock-run1-preit50.dat
if [ ! -s $newinp ]; then
 python $ATTRACTTOOLS/swapcombine.py $topstruc > $topstruc-combined
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $topstruc-combined --pdb proteinA.pdb proteinB.pdb $topstruc-gvm1
 python $ATTRACTTOOLS/fill-energies.py $topstruc-combined $topstruc-gvm1 > $topstruc-gvm2
 python $ATTRACTTOOLS/topcombined.py $topstruc-gvm2 $nstruc > $topstruc-topcombined
 python $ATTRACTTOOLS/monte.py $topstruc-topcombined $noclone > $newinp
fi
inp=$newinp

outp=dock-run1-it50.dat
sorted=dock-run1-it50-sorted.dat
topstruc=dock-run1-it50-topstruc.dat
if [ ! -s $outp-filtered ]; then
 python ../neo-attract-em.py $inp $dock0 $dock_stage4 --chunks $threads  --np $threads --output $outp
 python ../neo-attract-em.py $outp $dock0 $dock_stage4_min --chunks $threads  --np $threads --output $outp-min
 python $ATTRACTTOOLS/filter-energy.py $outp-min $energy_threshold > $outp-filtered
 
fi
if [ ! -s $topstruc ]; then
 python /home/cschindler/software/attract/tools/saxs/saxs-score-wrapper.py exp-processed.dat $outp-filtered --pdb proteinA.pdb proteinB.pdb  $outp-gvm1
 python $ATTRACTTOOLS/fill-energies.py $outp-filtered $outp-gvm1 > $outp-gvm2
 python $ATTRACTTOOLS/sort.py $outp-gvm2 > $sorted 
 $ATTRACTDIR/fix_receptor $sorted 2 | python $ATTRACTTOOLS/fill.py /dev/stdin $sorted > $sorted-fixre 
 $ATTRACTDIR/fix_receptor taboo.dat 2 | python $ATTRACTTOOLS/fill.py /dev/stdin taboo.dat > taboo.dat-fixre 
 python $ATTRACTDIR/dump_coordinates.py $sorted-fixre proteinB-CA.pdb structurescoor.npy -1 
 python $ATTRACTDIR/dump_coordinates.py taboo.dat-fixre proteinB-CA.pdb taboocoor.npy -1 
 python $ATTRACTTOOLS/saxs/savechi.py $sorted chiscore.npy 
 python $ATTRACTTOOLS/saxs/taboocluster.py 3.0 6.0 0 structurescoor.npy taboocoor.npy chiscore.npy taboochi.npy $sorted taboo.dat > $outp-filtered-tabu 
  
 python $ATTRACTTOOLS/sort.py $outp-filtered-tabu > $outp-filtered-tabu-sorted 
  $ATTRACTTOOLS/top $outp-filtered-tabu-sorted $ntop > $topstruc
fi

if [ ! -s result-run1.dat ]; then
 ln -s taboo.dat result-run1.dat
fi
