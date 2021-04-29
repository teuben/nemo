#
CC  = gcc
CXX = g++
FC  = gfortran

CFLAGS   = -g -O2
CXXFLAGS = -g -O2
FFLAGS   = -g -O2

TIME = /usr/bin/time
DAT  = $(NEMODAT)/nbabel/input256

BINS  = nbabel_c nbabel_cpp
BENCH = bench_c  bench_cpp

help:
	@echo "nbabel"

all:	$(BINS)

clean:  $(BINS)
	rm -f $(BINS)

bench:	$(BENCH)

nbabel_c:	nbabel.c
	$(CC) $(CFLAGS) -o nbabel_c nbabel.c -lm

nbabel_cpp:	main.cpp
	$(CXX) $(CXXFLAGS) -o nbabel_cpp main.cpp -lm



bench_c: nbabel_c
	$(TIME) ./nbabel_c < $(DAT)

bench_cpp: nbabel_cpp
	$(TIME) ./nbabel_cpp < $(DAT)
