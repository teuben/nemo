#
include $(NEMOLIB)/makedefs
#	some routines exist in single as well as double precision


#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory

BINFILES = cpu-bench_f cpu-bench_d

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`


clean:
	rm -f *.o *.a core $(BINFILES)

all:	$(BINFILES)

cpu-bench_f:
	$(CC) -DSINGLEPREC $(CFLAGS) -o cpu-bench_f  cpu-bench.c $(NEMO_LIBS)

cpu-bench_d:
	$(CC) -DDOUBLEPREC $(CFLAGS) -o cpu-bench_d  cpu-bench.c $(NEMO_LIBS)
