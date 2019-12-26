#	Makefile for NEMO/src/orbit/misc
include $(NEMOLIB)/makedefs

L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = dopri5.c dop853.c 
OBJFILES=  dopri5.o dop853.o
LOBJFILES= $L(dopri5.o) $L(dop853.o)
BINFILES = orbfour orbint orbintv orbplot otos perorb stoo orbsos orbstat orblist orboom
TESTFILES=  orbdim orblist

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

orbplot: orbplot.c
	$(CC) $(CFLAGS) -o orbplot orbplot.c $(NEMO_LIBS) $(YAPPLIB) $(EL) 

orboom: orboom.c
	$(CC) $(CFLAGS) $(HDF5_INC) -o orboom orboom.c $(NEMO_LIBS) $(HDF5_LIB) $(EL) 

