echo '**************************************************************'
echo 'cluster'
echo '**************************************************************'
$ATTRACTDIR/matrix-lrmsd out_$name.dat protein.pdb peptide/peptide.pdb --ens 2 peptide/peptides.list > out_$name.lrmsdlist

# parameters < RMSD cutoff in A> <minimum size of cluster>
$ATTRACTDIR/cluster_struc out_$name.lrmsdlist 5.0 1 > out_$name.cluster
python $ATTRACTDIR/cluster2dat.py out_$name.cluster out_$name.dat > out_$name-cluster.dat
cat out_$name-cluster.dat |awk '{print $4}' |sort > out_$name-cluster-centers
cat out_$name-cluster.dat |awk '{print $5}' |sort > out_$name-cluster-firsts
