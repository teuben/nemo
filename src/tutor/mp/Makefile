#

#   set defaults in case the makedefs does not
CC    = gcc
MPICC = mpicc

TIME = /usr/bin/time -f "%U %S %e %P"

include $(NEMOLIB)/makedefs

BINS = hello1 hello2 heated_plate_openmp poisson_mpi scaling scaling2

all:  $(BINS)

clean:
	rm -f $(BINS)

hello1:
	$(CC) -fopenmp -O3  -o hello1 hellomp.c  -lm

hello2:
	$(CC)          -O3  -o hello2 hellomp.c  -lm

bench: $(BINS)
	$(TIME) ./hello1
	$(TIME) ./hello2


heated_plate_openmp:
	$(CC) $(CFLAGS) -fopenmp -o heated_plate_openmp heated_plate_openmp.c $(NEMO_LIBS)

bench2: heated_plate_openmp
	$(TIME) ./heated_plate_openmp nprocs=4


# https://people.sc.fsu.edu/~jburkardt/c_src/poisson_mpi/poisson_mpi.c

poisson_mpi:
	$(MPICC) poisson_mpi.c -o poisson_mpi -lm

bench3:
	$(TIME) mpirun -np 4 ./poisson_mpi 512

# takes about 25" in single mode
bench4:
	$(TIME) ./sections iter=10000

# scaling from https://stackoverflow.com/questions/19780554/what-limits-scaling-in-this-simple-openmp-program

scaling: scaling.c
	$(CC) scaling.c -std=c99 -fopenmp -O3    -o scaling -lm

ONT=1 2 4 8 16
bench5: scaling
	@echo OMP_NUM_THREADS ONT=$(ONT)
	-@for i in $(ONT); do\
	   (echo $$i;OMP_NUM_THREADS=$$i $(TIME) ./scaling); done


quad_openmp: quad_openmp.c
	gcc -fopenmp quad_openmp.c -o quad_openmp

bench6:  quad_openmp
	@echo OMP_NUM_THREADS ONT=$(ONT)
	-@for i in $(ONT); do\
	   (echo $$i;OMP_NUM_THREADS=$$i $(TIME) ./quad_openmp); done


bench7: python_openmp.py
	@echo OMP_NUM_THREADS ONT=$(ONT)
	-@for i in $(ONT); do\
	   (echo $$i;OMP_NUM_THREADS=$$i $(TIME) ipython python_openmp.py); done

bench8: scaling2
	@echo "# OMP_NUM_THREADS ONT=$(ONT)"
	-@for i in $(ONT); do\
	   (echo -n "$$i ";$(TIME) ./scaling2 umax=20000 np=$$i iter=20*$$i debug=-1); done


bench9:	scaling2
	@echo Long integration, see laptop performance drop, see core swaps
	$(TIME) ./scaling2  iter=500 > bench9.log 2>&1
	grep cputime bench9.log | sed s/###// | tabplot - 6 8 line=1,1 ycoord=0  ymin=0
	@echo "Jansky 0.077  k2 0.090 0.098"
