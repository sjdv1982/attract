#!/usr/bin/bash
rm resultk.dat
gawk -f one.awk va1=1 start.pdb > on.pdb
pdb4amber -i on.pdb -y -d -p -o one.pdb
gawk -f two.awk va1=1 start.pdb > tw.pdb
pdb4amber -i tw.pdb -y -d -p -o two.pdb
tleap -f one_leap.in
tleap -f two_leap.in
natom=`gawk -f atom.awk one.crd`
echo $natom
c="0"
while [ $c -lt 1000 ]
do
c=$[$c+1]
echo "refinement of structure " $c
gawk -f get2.awk va1=$c start.pdb > mu.pdb
python ~/scripts/rename_pdb_for_chimera.py 2 mu.pdb
pdb4amber -i mu-collect.pdb -y -d -p -o mur.pdb
gawk -f ter.awk mur.pdb > murks2.pdb
gawk -f ter2.awk murks2.pdb > murks.pdb
tleap -f  leap.in
echo "Start minimization"
timeout 120 bash mini.sh
if [ $? -eq 124 ]; then
    echo "Minimization failed, copying existing files"
    if [ -e exp1.crd ]; then
         cp exp1.crd exp2.crd
    elif [ -e exp0.crd ]; then
        cp exp0.crd exp2.crd
    else
       cp exp.crd exp2.crd
    # Timeout occurred
    fi
fi

ambpdb -p exp.top < exp2.crd > exp2.pdb
./inter exp2.pdb 1 > dist.in
cp exp2.pdb 'bs'$c'.pdb'
cp exp2.crd 'bs'$c'.crd'
pmemd.cuda -O -i md2.in -c exp2.crd -p exp.top -r exp3.crd -o exp3.out
pmemd.cuda -O -i md3.in -c exp3.crd -p exp.top -r exp4.crd -o exp4.out
pmemd.cuda -O -i minc.in -c exp4.crd -p exp.top -r exp5.crd -o exp5.out
pmemd.cuda -O -i minc.in -c exp5.crd -p exp.top -r exp6.crd -o exp6.out
./onetwo exp6.crd $natom one1.crd two1.crd
mpiexec -np 4 sander.MPI -O -i mingb.in -c exp6.crd -p exp.top -r exp7.crd -o exp7.out
mpiexec -np 4 sander.MPI -O -i mingb.in -c one1.crd -p one.top -r one2.crd -o one2.out
mpiexec -np 4 sander.MPI -O -i mingb.in -c two1.crd -p two.top -r two2.crd -o two2.out
va1=`gawk -f ener.awk exp7.out`
va2=`gawk -f ener.awk one2.out`
va3=`gawk -f ener.awk two2.out`
va4=`gawk -f eave.awk exp4.out`
echo $c $va1 $va2 $va3 $va4 >> resultk.dat
mv -f exp6.crd 'bsf'$c'.crd'
ambpdb -p exp.top < 'bsf'$c'.crd' > 'bsf'$c'.pdb'
done
