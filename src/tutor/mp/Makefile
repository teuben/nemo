#

CC = gcc

TIME = /usr/bin/time

include $(NEMOLIB)/makedefs

BINS = hello1 hello2 heated_plate_openmp poisson_mpi scaling

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
	mpicc poisson_mpi.c -o poisson_mpi -lm

bench3:
	$(TIME) mpirun -np 4 ./poisson_mpi 512

# takes about 25" in single mode
bench4:
	$(TIME) ./sections iter=10000

# scaling from https://stackoverflow.com/questions/19780554/what-limits-scaling-in-this-simple-openmp-program

scaling: scaling.c
	$(CC) scaling.c -std=c99 -fopenmp -O3    -o scaling

ONT=1 2 4 8
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


