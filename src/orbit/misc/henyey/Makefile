#	Testing FORTRAN interface in HENYEY

include $(NEMOLIB)/makedefs

TIME = /usr/bin/time

#  P is used by 'forcon', to convert code from REAL*4 <---> REAL*8
#  P must be 4 or 8 when called (e.g. "make real P=8" makes a REAL*8 version)
P=8
#  B is the bar number (B=2 would give bar2 to be added to henyey)
B=6
#  Compiler flags
CFLAGS = -g -DREAL$(P)
FFLAGS = -g -std=legacy  -fallow-argument-mismatch

help:
	@echo "Henyey targets"
	@echo "   clean "
	@echo "   src    into  SRCDIR=$(SRCDIR)"
	@echo "   henyey B=$(P)"
	@echo "	  real P=$(P)"
	@echo "   tar"

SRCS = Makefile README henyey.f henyey.h maxorb.h model.h forcon

BINS = henyey bar bar2 bar6 

SRCDIR = $(NEMO)/src/orbit/misc/henyey

install:
	@echo "Not part of install yet"

clean:
	rm -f *.o *.a core $(BINS) fort.20 \
		test_main.c henyey_main.c bar_main.c

bin:	$(BINS)

src:	.src
.src:	$(SRCS)
	cp $? $(SRCDIR)
	@touch .src

henyey:	henyey.o bar$(B).o henyey_main.o
	$(FC) $(FFLAGS) -o henyey henyey.o henyey_main.o bar$(B).o \
	$(NEMO_LIBS) $(LOCAL_LIB) $(FORLIBS) $(EL)

henyey_main.o: 
	$(CC) $(CFLAGS) -I$(NEMOINC) -I$(NEMOLIB) -c henyey_main.c

henyey_main.c: henyey.f
	ftoc henyey.f > henyey_main.c

henyey.o:  henyey.f
	ftoc henyey.f > henyey_main.c
	$(FC) $(FFLAGS) -c henyey.f


bar:	bar.o
	$(FC) $(FFLAGS) -o bar bar.o bar_main.o \
	  -L$(NEMOLIB) -lnemo

bar.o:  bar.f
	ftoc bar.f > bar_main.c
	$(CC) $(CFLAGS) -c bar_main.c
	$(FC) $(FFLAGS) -c bar.f

bar6:   bar6.c
	$(CC) -DTESTBED $(CFLAGS) -o bar6 bar6.c -lm

bar2:   bar2.c
	$(CC) -DTESTBED $(CFLAGS) -o bar2 bar2.c -lm

real:
	@echo "Making a real*$(P) version of the code (P=$(P))"
	./forcon -$(P) henyey.f *.h
	@echo "$(P)" > P

# no rotation, launch from X axis
bench1:
	$(TIME) ./henyey 0.6 0 0 1.9 1.0 1 0.0 0

bench1m:
	$(TIME) ./henyey 0.6 0 0 1.9 1.0 1 0.0 1000 0.001 > tab1
	grep ^'  1  0.00' tab1 > tab1x
	grep ^'201  '     tab1 > tab1y

# no rotation, launch from Y axis
bench2:
	./henyey 0 0.6 1.1 0 1.1 2 0.0 0

# some rotation
bench3:
	./henyey 0.6 0 0 1 1.3 1 0.2 0

# as bench1, but using perorb
bench4:
	perorb phase=0.6,0,0,0,1.85,0 dir=x potname=bar83 potpars=0,1,1,0.2,1 debug=1

bench4m:
	$(TIME) perorb phase=0.6,0,0,0,1.85,0 dir=x potname=bar83 potpars=0,1,1,0.2,1  step=0.001 norbit=1000 > tab4


