#	Makefile for NEMO/src/nbody/trans
include $(NEMOLIB)/makedefs

#CFLAGS = -O
#LDFLAGS = -O
#RANLIB = ranlib

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
SRCFILES= 
OBJFILES=
LOBJFILES=
BINFILES = rungalaxy
TESTFILES= 

help:
	@echo NEMO 

clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	update $(BINFILES) $(TESTFILES)
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
#install:   .install_inc .install_lib 
install:	$(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)
	@touch .install_lib

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)
	@touch .install_inc


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

# Targets used to assemble object library.

bin: $(BINFILES)
tests: $(TESTFILES)
