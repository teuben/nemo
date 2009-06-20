#	Makefile for NEMO/src/nbody/image
include $(NEMOLIB)/makedefs

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = snapccd.c snapgrid.c snapmap.c snapslit.c
OBJFILES=  
LOBJFILES= 
BINFILES = snapccd snapgrid snapmap snapslit
TESTFILES= 

help:
	@echo NEMO V2.2 NEMO/src/nbody/image

clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

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

#install:   .install_inc .install_lib 
install:

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

bin: $(BINFILES)
tests: $(TESTFILES)

# special

snapslit:  snapslit.c
	$(CC) $(CFLAGS) -o snapslit snapslit.c $(NEMO_LIBS) $(YAPPLIB) -lm

