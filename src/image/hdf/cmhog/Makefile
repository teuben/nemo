#	Makefile for NEMO/src/kernel/hdf
include $(NEMOLIB)/makedefs
LOCAL_INC = $(HDF_INC)
LOCAL_LIB = $(HDF_LIB)

#	sub-directories to be visited, in the order given here
#DIRS = 

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = hdfgrid.c
OBJFILES=  
LOBJFILES= 
BINFILES = 
TESTFILES= 

ifneq ($(strip $(HDF_LIB)),)
  BINFILES += hdfgrid
endif

help:
	@echo NEMO/src/image/hdf/cmog: Support files for the cmhog hydro project
	@echo OBJFILES=$(OBJFILES)
	@echo bin: BINFILES=$(BINFILES)
	@echo tests: TESTFILES=$(TESTFILES)
	@echo HDF_INC=$(HDF_INC)
	@echo HDF_LIB=$(HDF_LIB)
	@echo 'on e.g. linux you also need CLFLAGS="-DPC -DPC386 -Dunix"'
clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

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

install:   

.install_lib: $(OBJFILES) 
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)
	@touch .install_inc

.install_bin: $(BINFILES)
	mv $? $(NEMOBIN)

# Targets used by Nemo to update libraries and bins.

nemo_lib: $(LIBFILE)

nemo_sublibs:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_lib); done

nemo_bin: nemo_subbins $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	rm -f *.o

nemo_subbins:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_bin); done

nemo_src:
	-@for i in $(SUBDIRS); do \
	(cd $$i; $(MAKE) nemo_src); done
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

#	Special rules
bin: $(BINFILES)
tests: $(TESTFILES)

#	only do this at Maryland
tar:
	@echo Not activated yet
#	(cd $(NEMO); tar cvf /home/ftp/progs/nemo/cmhog.tar usr/pjt/cmhog)
