#

help:
	@echo "Heggie & Hut 2003 : The Gravitational Million-Body Problem"

all:	nbody_hh

nbody_hh: nbody_hh.c
	$(CC) $(CFLAGS) -o nbody_hh nbody_hh.c -lm

bench:
	@echo Benchmark

clean:
	rm *.o nbody_hh
