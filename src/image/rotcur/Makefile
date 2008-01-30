#	Makefile for NEMO/src/kernel/io
include $(NEMOLIB)/makedefs

## **************** Attempt for template ************************

#	Compile flags (the default in this dir is single precision images
PRECFLAG = 
# CFLAGS = -O $(PRECFLAG)
# LDFLAGS = -O $(PRECFLAG)
# RANLIB = ranlib

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
SRCFILES = 
OBJFILES=  rotcurs.o
LOBJFILES= $L(rotcurs.o)
BINFILES = ccdvel rotcur rotcurshape rotcurtab pvtrace velcube velfit
TESTFILES= 

help:
	@echo NEMO/src/image/rotcur

clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

#	update or install a library
update:  .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

install: .install_lib

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



# special


ccdplot:  ccdplot.c
	$(CC) $(CFLAGS) -o ccdplot ccdplot.c \
		$(L) $(YAPPLIB) $(EL) -lm
