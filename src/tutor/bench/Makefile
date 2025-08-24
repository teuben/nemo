#
include $(NEMOLIB)/makedefs
#	some routines exist in single as well as double precision

# gnu time
TIME = /usr/bin/time -f '%Uuser %Ssystem %Eelapsed %PCPU'


#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory

BINFILES = cpu-bench_f cpu-bench_d big big8 offt offt8 flops


.PHONY : help
## help:      This Help
help : Makefile
	@sed -n 's/^##//p' $<


clean:
	rm -f *.o *.a core $(BINFILES)

all:	$(BINFILES)

## flops:    old MFLOP estimate
flops:	flops.c
	$(CC) -DPOSIX1 $(CFLAGS) -o flops flops.c -lm

## cpu-bench_f:
cpu-bench_f:
	$(CC) -DSINGLEPREC $(CFLAGS) -o cpu-bench_f  cpu-bench.c $(NEMO_LIBS)

## cpu-bench_d:
cpu-bench_d:
	$(CC) -DDOUBLEPREC $(CFLAGS) -o cpu-bench_d  cpu-bench.c $(NEMO_LIBS)

## bench:  runs cpu-bench_{f,d} and flops
bench:	$(BINFILES)
	time ./cpu-bench_f > cpu-bench_f.log
	time ./cpu-bench_d > cpu-bench_d.log
	time ./flops 

## bench8:
bench8:	big big8 offt offt64
	@./big | wc
	@./big8 | wc
	offt
	offt8


NP = 1 2 4 6 8 10 12 14 16 18 20 22 24
## bench10:   sdmath for bunch of NP's via NP=
bench10:
	rm -f bench10.1.log bench10.2.log
	for np in $(NP); do\
	   echo -n ":: $$np :: ";\
           echo "$$np " >> bench10.1.log;\
	   $(TIME) sdmath iter=100 np=$$np 2>>bench10.2.log;\
	done
	@paste bench10.1.log bench10.2.log > bench10.log
	@echo "See bench10.log for results"
	@../../scripts/python/amdahl.py bench10.log

#
big:	big.f
	$(FC) -o big big.f

big8:	big8.f
	$(FC) -o big8 big8.f

offt:	offt.c
	$(CC) -o offt offt.c

offt8:	offt.c
	$(CC) -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -o offt8 offt.c


array_add_omp:  array_add.c
	$(CC) $(CFLAGS)  -fopenmp array_add.c -o array_add_omp $(NEMO_LIBS)
#	$(CC) $(CFLAGS) -std=c99 -fopenmp array_add.c -o array_add_omp $(NEMO_LIBS)
