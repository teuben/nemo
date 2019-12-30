#
include $(NEMOLIB)/makedefs
#


L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = 
OBJFILES=  
LOBJFILES= 
BINFILES = runlpvi
TESTFILES=  

help:
	@echo NEMO V2.2 NEMO/src/orbit/misc
	@echo BINFILES=$(BINFILES)

clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)
bin: $(BINFILES)
#	update or install a library
update:  .update_inc .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

# install:   .install_inc .install_lib 
install: .install_lib

.install_lib: $(OBJFILES) 
	ar ruv $(L) $?
	$(RANLIB) $(L)

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

nemo_subbins:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_bin); done

nemo_src:
	-@for i in $(SUBDIRS); do \
	(cd $$i; $(MAKE) nemo_src); done
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

#	Special rules

tests: $(TESTFILES)

#	

orboom: orboom.c
	$(CC) $(CFLAGS) $(HDF5_INC) -o orboom orboom.c $(NEMO_LIBS) $(HDF5_LIB) $(EL) 

runlpvi: runlpvi.c
	$(CC) $(CFLAGS) -o runlpvi runlpvi.c $(NEMO_LIBS) $(EL) 




run1:
	runlpvi run1 posvel=0.2,0.34,0.5,0.1235 prefix=log dt=0.05 tstop=25000 orbit=1 ci=1::7 nstep=20
	tabplot run1/log.orb 2 3 -0.5 0.5 -0.5 0.5 yapp=1/xs
	mkorbit run1/orb.in 0.2 0.34 0 0.5 0.1235 0 potname=log potpars=0,2,0,0.7
	orbint run1/orb.in run1/orb.out  50000 0.05 ndiag=10000
	orbplot run1/orb.out xrange=-0.5:0.5 yrange=-0.5:0.5 yapp=2/xs

orb:
	rm run1/orb.*
	mkorbit run1/orb.in 0.2 0.34 0 0.5 0.1235 0 potname=log potpars=0,2,0,0.7
	orbint run1/orb.in run1/orb.out  500 0.05 ndiag=100
	orbplot run1/orb.out xrange=-0.5:0.5 yrange=-0.5:0.5 yapp=2/xs


