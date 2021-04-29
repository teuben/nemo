#
CC = gcc
CFLAGS = -g -O2

TIME = /usr/bin/time

help:
	@echo "Heggie & Hut 2003 : The Gravitational Million-Body Problem"

all:	nbody_hh

nbody_hh: nbody_hh.c
	$(CC) $(CFLAGS) -o nbody_hh nbody_hh.c -lm

bench:
	@echo Benchmark
	$(TIME) ./nbody_hh 1024 10 1.0

clean:
	rm -f *.o nbody_hh
