#	Testing FORTRAN interface in HENYEY

#  P is used by 'forcon', to convert code from REAL*4 <---> REAL*8
#  P must be 4 or 8 when called (e.g. "make real P=8" makes a REAL*8 version)
P=4
#  B is the bar number (B=2 would give bar2 to be added to henyey)
B=6
#  Compiler flags
CFLAGS = -g -DREAL$(P)
FFLAGS = -g

help:
	@echo Henyey targets
	@echo "   clean "
	@echo "   src    into  SRCDIR=$(SRCDIR)"
	@echo "   henyey B=$(P)"
	@echo "	  real P=$(P)
	@echo "   tar"
	

SRCS = Makefile README henyey.f henyey.h maxorb.h model.h \
	forcon

BINS = henyey bar bar2 bar6 

SRCDIR = $(NEMO)/src/orbit/misc/henyey

install:
	@echo Not part of install yet

clean:
	rm -f *.o *.a core $(BINS) fort.20 \
 		test_main.c henyey_main.c bar_main.c

bin:	$(BINS)

src:	.src
.src:	$(SRCS)
	cp $? $(SRCDIR)
	@touch .src

henyey:	henyey.o bar$(B).o
	$(FC) $(FFLAGS) -o henyey henyey.o henyey_main.o bar$(B).o \
	  -L$(NEMOLIB) -lnemo

henyey.o:  henyey.f
	ftoc henyey.f > henyey_main.c
	$(CC) $(CFLAGS) -c henyey_main.c
	$(FC) $(FFLAGS) -c henyey.f


bar:	bar.o
	$(FC) $(FFLAGS) -o bar bar.o bar_main.o \
	  -L$(NEMOLIB) -lnemo

bar.o:  bar.f
	ftoc bar.f > bar_main.c
	$(CC) $(CFLAGS) -c bar_main.c
	$(FC) $(FFLAGS) -c bar.f

bar6:
	$(CC) -DTESTBED $(CFLAGS) -o bar6 bar6.c -lm

bar2:
	$(CC) -DTESTBED $(CFLAGS) -o bar2 bar2.c -lm

tar:
	tar cvf henyey.tar $(SRCS)

real:
	@echo "Making a real*$(P) version of the code (P=$(P))"
	forcon -$(P) henyey.f *.h
	@echo "$(P)" > P
