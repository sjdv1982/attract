#!/usr/bin/bash
rm exp1.crd
rm exp2.crd
rm exp0.crd
mpiexec -np 4 sander.MPI -O -i mi0.in -c exp.crd -p exp.top -r exp0.crd -o exp0.out
mpiexec -np 4 sander.MPI -O -i mi.in -c exp0.crd -p exp.top -r exp1.crd -o exp1.out
mpiexec -np 4 sander.MPI -O -i mi2.in -c exp1.crd -p exp.top -r exp2.crd -o exp2.out