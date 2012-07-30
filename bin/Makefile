OBJ= read_hm.o read_dof.o read_parameters.o read_pdb.o write_pdb.o print_struc.o \
  cartstate.o euler2rotmat.o rotate.o deform.o minfor.o monte.o mc11a.o GGUBS.o\
  parse_options.o select.o pairgen.o pairenergy.o energy.o globalenergy.o \
  matmult.o matinv.o prox.o reset_forces.o nonbon8.o ligmin.o rota.o trans.o \
  em/situs/lib_vio.o em/situs/lib_pio.o em/situs/lib_err.o \
  parse_restraints.o restrain.o read_ens.o read_pdb2.cpp\
  axisrot.o crand.o sym.o enstrans.o ene_morph.cpp \
  euler2torquemat.o memcpy.o disre.o moderest.o
    
FLAG= -g -p -O3 -fno-automatic -ffast-math -fcray-pointer
#FLAG= -g -p -O0 -fno-automatic -ffast-math -fcray-pointer
H=state.h max.h grid.h nonbon.h makegrid.h prox.h
FF=gfortran

CXX        = g++
CC        = gcc
CFLAGS = -g -p -O3 -Wall -ffast-math 
#CFLAGS = -g -p -O0 -fno-inline -Wall -ffast-math 

attract: attract.o em/emenergy.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@

attract-torque: attract.o em/emenergy.o grid-torque.o ministate-torque.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@

attract-em-gradcheck: attract.o em/emenergy-gradcheck.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@

attract-em-precompcheck: attract.o em/emenergy-precompcheck.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@

attract-em-contour: attract.o em/emenergy-contour.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@

collect: collect.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o 
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

emcoor: emcoor.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o em/lib_em.o em/lib_corr.o em/situs/lib_vio.o em/situs/lib_pio.o em/situs/lib_err.o
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

gvm: gvm.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o em/lib_em.o em/lib_corr.o em/situs/lib_vio.o em/situs/lib_pio.o em/situs/lib_err.o
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

matrix-lrmsd: matrix-lrmsd.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o 
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

lrmsd: lrmsd.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o 
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

rmsd: rmsd.o read_parameters.o read_pdb.o read_pdb2.o read_dof.o \
  cartstate.o euler2rotmat.o deform.o rotate.o read_hm.o memcpy.o write_pdb.o \
  read_ens.o
	$(CXX) $(CFLAGS)  $^ -lgfortran -o $@

calc_interior: calc_interior.o em/situs/lib_vio.o em/situs/lib_pio.o em/situs/lib_err.o
	$(CXX) $(CFLAGS) $^  -lgfortran -o $@

make-grid: make-grid.o em/emenergy.o grid_calculate.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS) $^  -lgfortran -lrt -o $@

make-grid-omp: make-grid.o em/emenergy.o grid_calculate-omp.o grid.o ministate.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^  -fopenmp -lgfortran -lrt -o $@

make-grid-torque: make-grid-torque.o em/emenergy.o grid_calculate-torque.o grid-torque.o ministate-torque.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^  -lgfortran -lrt -o $@

make-grid-torque-omp: make-grid-torque.o em/emenergy.o grid_calculate-torque-omp.o grid-torque.o ministate-torque.o $(OBJ) 
	$(CXX) $(CFLAGS)  $^  -fopenmp -lgfortran -lrt -o $@

shm-grid: shm-grid.o grid_calculate.o grid.o ministate.o em/emenergy.o $(OBJ)
	$(CXX) $(CFLAGS) $^  -lgfortran -lrt -o $@

shm-grid-torque: shm-grid-torque.o grid_calculate-torque.o grid-torque.o ministate-torque.o em/emenergy.o $(OBJ)
	$(CXX) $(CFLAGS)  $^ -lgfortran -lrt -o $@ 

shm-clean: shm-clean.o
	$(CXX) $(CFLAGS)  $^ -lrt -o $@
	
systsearch: systsearch.o
	$(FF) $(FLAG) $< -o $@
	
grid.o: grid.cpp $(H) 
	$(CXX) $(CFLAGS) --param max-inline-insns-single=2000 -c $< -o $@ 

grid-torque.o: grid.cpp $(H) 
	$(CXX) $(CFLAGS) --param max-inline-insns-single=2000 -DTORQUEGRID -c $< -o $@ 

map2map:	em/situs/map2map.o em/situs/lib_err.o em/situs/lib_vio.o \
		em/situs/lib_vwk.o em/situs/lib_vec.o em/situs/lib_std.o 
	$(CC) $(CFLAGS) $^ -o $@ -lm

fix_receptor:   fix_receptor.cpp read_dof.o euler2rotmat.o matinv.o matmult.o \
                print_struc.o max.h
	$(CXX) $(CFLAGS) $^ -o $@ -lm

deredundant: deredundant.cpp read_dof.o euler2rotmat.o matmult.o \
                print_struc.o max.h
	$(CXX) $(CFLAGS) $^ -o $@ -lm

omp_threadnum: omp_threadnum.cpp
	$(CXX) $(CFLAGS)  $^  -fopenmp -o $@

%-torque-omp.o: %.cpp $(H)
	$(CXX) $(CFLAGS) -fopenmp -DTORQUEGRID -c $< -o $@ 

%-torque.o: %.cpp $(H)
	$(CXX) $(CFLAGS) -DTORQUEGRID -c $< -o $@ 

%-omp.o: %.cpp $(H)
	$(CXX) $(CFLAGS) -fopenmp -c $< -o $@ 

em/emenergy-contour.o:  em/emenergy.cpp
	$(CXX) $(CFLAGS) -DCONTOUR -c $< -o $@

em/emenergy-gradcheck.o:  em/emenergy.cpp
	$(CXX) $(CFLAGS) -DGRADCHECK -c $< -o $@

em/emenergy-precompcheck.o:  em/emenergy.cpp
	$(CXX) $(CFLAGS) -DPRECOMPCHECK -c $< -o $@

collectlib.so: collectlib-pic.o read_parameters-pic.o read_pdb-pic.o \
  read_pdb2-pic.o read_dof-pic.o cartstate-pic.o euler2rotmat-pic.o \
  deform-pic.o rotate-pic.o read_hm-pic.o memcpy-pic.o write_pdb-pic.o \
  read_ens-pic.o
	$(CXX) $(CFLAGS) -Wall -shared -fPIC $^ -lgfortran -lm -o $@ 

collectlib: collectlib.so
	
%.o: %.cpp $(H)
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.c $(H)
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.f
	$(FF) $(FLAG) -c $< -o $@

%-pic.o: %.cpp $(H)
	$(CXX) $(CFLAGS) -c -fPIC $< -o $@

%-pic.o: %.f
	$(FF) $(FLAG) -c -fPIC $< -o $@

clean: 
	rm -f *.o
	rm -f em/*.o 
	rm -f em/situs/*.o 
	rm -f attract

cluster_struc: cluster_struc.o  
	$(CXX) $(CFLAGS)  $^ -o $@
	
rest: read_pdb2.o
	$(FF) $(FLAG) reduce.f -o reduce
	$(FF) $(FLAG) translate.f -o translate 
	$(FF) $(FLAG) rotam.f -o rotam 
	$(FF) $(FLAG) modes.f -o modes 
	$(FF) $(FLAG) modesca.f -o modesca 
	$(FF) $(FLAG) compare.f -o compare 
	$(FF) $(FLAG) viewe.f -o viewe 
	$(FF) $(FLAG) rmsca.f -o rmsca 
	$(CXX) $(CFLAGS) center.cpp read_pdb2.o -o center
	$(CXX) $(CFLAGS) transl.cpp read_pdb2.o -o trans
	$(CXX) $(CFLAGS) deflate.cpp read_pdb2.o -o deflate
        
all: attract attract-torque collect lrmsd rmsd calc_interior make-grid \
 make-grid-omp make-grid-torque make-grid-torque-omp shm-grid \
 shm-grid-torque shm-clean systsearch map2map fix_receptor \
 omp_threadnum deredundant matrix-lrmsd rest 

