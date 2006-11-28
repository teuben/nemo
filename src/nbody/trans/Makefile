#	Makefile for NEMO/src/nbody/trans
include $(NEMOLIB)/makedefs

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory

MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES= 
OBJFILES=
LOBJFILES=
BINFILES = snapcenter snapsort snapstack snaptrim snaprotate snapdens \
	snapscale unbind snapmask snapadd snaprect \
        snapsphere snapmass snapspin snaptrans snapvirial \
        snapcopy snapinert snapmerge snapshift snapsplit \
        snapmerge_a_sp snapmerge_a_dp

TESTFILES =

help:
	@echo NEMO V2.2 NEMO/src/nbody/trans

clean:	cleanbench
	-rm -f *.o *.a core $(BINFILES) $(TESTFILES)
	-(cd hackdens; make clean)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)
bin: $(BINFILES)
test: $(TESTFILES)
#	update the library: direct comparison with modules inside L
update:  .update_inc .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

#	install the library: can be (re)done from scratch by removing .o files

install:

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)
	@touch .install_lib

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)
	@touch .install_inc

# Targets used by Nemo to update libraries and bins.

nemo_lib: $(LIBFILE)

nemo_sublibs:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_lib); done

nemo_bin: nemo_subbins $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES) 
	mv $(BINFILES) $(NEMOBIN)
	(cd hackdens; $(MAKE) bins)

nemo_subbins:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_bin); done

nemo_src:
	-@for i in $(SUBDIRS); do \
	(cd $$i; $(MAKE) nemo_src); done
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

# benchmark
ITER= 10
NP6 = 1000000

p6:
	mkplummer p6 $(NP6) massname='n(m)' massrange=1,2


bench1: snapbench p6
	@echo "using NP6=$(NP6) ITER=$(ITER), use 'make cleanbench' if you use a new NP6"
	time ./snapbench p6 'mass=3.1415'   bodytrans=f iter=$(ITER)
	time ./snapbench p6 'mass=3.1415*m' bodytrans=t iter=$(ITER)

cleanbench:
	rm -f p6

# Targets used to assemble object library.

bin: $(BINFILES)
tests: $(TESTFILES)

hackdens: .hackdens

.hackdens:
	(cd hackdens ; $(MAKE) hackdens)

# targets for snapmerge_a_(s|d)_p
snapmerge_a : snapmerge_a_sp snapmerge_a_dp
	@/bin/cp snapmerge_a_sp snapmerge_a_dp $(NEMOBIN)

snapmerge_a_sp : snapmerge_a_sp.o 
	$(CC) -o $@ snapmerge_a_sp.o  -L$(NEMOLIB) -lnemo -lm

snapmerge_a_sp.o : snapmerge_a.c
	$(CC) $(CFLAGS) -DSINGLEPREC -c snapmerge_a.c -o snapmerge_a_sp.o

snapmerge_a_dp : snapmerge_a_dp.o 
	$(CC) -o $@ snapmerge_a_dp.o  -L$(NEMOLIB) -lnemo -lm

snapmerge_a_dp.o : snapmerge_a.c
	$(CC) $(CFLAGS) -c snapmerge_a.c -o snapmerge_a_dp.o
#
