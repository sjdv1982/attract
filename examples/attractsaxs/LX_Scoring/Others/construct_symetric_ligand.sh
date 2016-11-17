#ln -s *A-unbound.pdb refeA.pdb
#ln -s *B-unbound.pdb refeB.pdb

symetric=$1 # "A" or "B" for reference A or B
symindex=$2 # There might be several symetries

ligand=A
if [ $symetric = A ];then ligand=B;fi

echo "##### make symmetric ligand#####"
center=`$ATTRACTDIR/center refe$symetric.pdb |& grep Center| awk '{print $2, $3, $4}'`
echo "#pivot 1 " $center > ../structure-single-sym.dat
echo "#centered receptor: false" >> ../structure-single-sym.dat
echo "#centered ligands: false" >> ../structure-single-sym.dat
echo "#1" >> ../structure-single-sym.dat

rot=`python $ATTRACTTOOLS/euler.py refe$symetric.pdb refe$symetric-sym$symindex.pdb`
echo $rot >> ../structure-single-sym.dat

$ATTRACTDIR/collect ../structure-single-sym.dat refe$ligand.pdb --single > refe$ligand-to-${symetric}sym$2.pdb

#echo '##### shift ligand ######'
#centerb=`$ATTRACTDIR/center refeA.pdb |& grep Center| awk '{print $2,  $3, $4}'`
#centera=`$ATTRACTDIR/center refeA*sym$2.pdb |& grep Center| awk '{print $2, $3, $4}'`
#echo "#pivot 1 " $centera > ../structure-single-sym1.dat
#echo "#centered receptor: false" >> ../structure-single-sym1.dat
#echo "#centered ligands: false" >> ../structure-single-sym1.dat
#echo "#1" >> ../structure-single-sym1.dat
#b=(`echo $centera | tr " " "\n"`)
#a=(`echo $centerb | tr " " "\n"`)
#x=`echo "${a[0]} ${b[0]}" | awk '{print $1-$2}'`
#y=`echo "${a[1]} ${b[1]}" | awk '{print $1-$2}'`
#z=`echo "${a[2]} ${b[2]}" | awk '{print $1-$2}'`
#echo "0.0 0.0 0.0 " $x $y $z >> ../structure-single-sym1.dat

#$ATTRACTDIR/collect ../structure-single-sym1.dat refeB*sym$2.pdb --single > refe$ligand-to-${symetric}sym$2-centered.pdb


echo "##### align structures #####"

center2=`$ATTRACTDIR/center refeA.pdb |& grep Center| awk '{print $2, $3, $4}'`
echo "#pivot 1 " $center2 > ../align-sym.dat
echo "#centered receptor: false" >> ../align-sym.dat
echo "#centered ligands: false" >> ../align-sym.dat
echo "#1" >> ../align-sym.dat
rot=`python $ATTRACTTOOLS/euler.py refeA.pdb *A-refe.pdb`
echo $rot >> ../align-sym.dat

$ATTRACTDIR/collect ../align-sym.dat refeB*sym$2.pdb --single > refeB-align-sym$2.pdb
$ATTRACTDIR/collect ../align-sym.dat refeA.pdb --single > refeA-align-sym$2.pdb
