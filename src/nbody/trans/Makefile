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
        snapcopy snapinert snapmerge snapshift snapsplit

TESTFILES =

help:
	@echo NEMO V2.2 NEMO/src/nbody/trans

clean:
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

p6:
	mkplummer p6 1000000 massname='n(m)' massrange=1,2

bench1: snapbench p6
	time ./snapbench p6 'mass=3.1415'   bodytrans=f iter=10
	time ./snapbench p6 'mass=3.1415*m' bodytrans=f iter=10 


# Targets used to assemble object library.

bin: $(BINFILES)
tests: $(TESTFILES)

hackdens: .hackdens

.hackdens:
	(cd hackdens ; $(MAKE) hackdens)


