#PBS -l nodes=1:ppn=12
#PBS -l walltime=6:30:00
#PBS -j oe
cd $PBS_O_WORKDIR

ATTRACTDIR=/home/sjoerd/data/work/attract

$ATTRACTDIR/bin/shm-clean

python $ATTRACTDIR/tools/randsearch.py 1 10000 > start.dat

$ATTRACTDIR/bin/calc_no_interior partner-1.pdb partner-1-interior.vol

$ATTRACTDIR/bin/reduce partner-1.pdb > /dev/null
echo 1 > partners.pdb
grep ATOM partner-1r.pdb >> partners.pdb
echo TER >> partners.pdb
echo END >> partners.pdb

$ATTRACTDIR/bin/make-grid partner-1r.pdb partner-1-interior.vol $ATTRACTDIR/parmu.par 10.000000 12.000000 partner-1.grid
$ATTRACTDIR/bin/shm-grid partner-1.grid partner-1.gridheader


python $ATTRACTDIR/protocols/attract.py start.dat $ATTRACTDIR/parmu.par partners.pdb --grid 1 partner-1.gridheader --gravity 1 --rstk 0.200000 --jobsize 1000 --np 12 --output dock.dat

$ATTRACTDIR/bin/shm-clean


