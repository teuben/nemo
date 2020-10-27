#

TIME = /usr/bin/time

include $(NEMOLIB)/makedefs

BINS = hello1 hello2 heated_plate_openmp poisson_mpi scaling

all:  $(BINS)

clean:
	rm -f $(BINS)

hello1:
	gcc -fopenmp -O3  -o hello1 hellomp.c  -lm

hello2:
	gcc          -O3  -o hello2 hellomp.c  -lm

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

scaling:
	gcc scaling.c -fopenmp -O3    -o scaling

bench5: scaling
	OMP_NUM_THREADS=1 $(TIME) ./scaling
	OMP_NUM_THREADS=2 $(TIME) ./scaling
	OMP_NUM_THREADS=4 $(TIME) ./scaling
	OMP_NUM_THREADS=8 $(TIME) ./scaling

