#	Makefile for NEMO/src/kernel/tab/fit
#       .so functions in here are currently only used by tabnllsqfit
include $(NEMOLIB)/makedefs

#	Compile flags  $(NFLAGS)
DIM=THREEDIM
NFLAGS=-D$(DIM)
LOCAL_INC = $(NFLAGS) -Iinc -fpic

#	remote dir where object files should live
R = $(NEMOOBJ)/fit

#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# testing a single fit
TEST = 

# Files etc. in this directory

MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 

FIT_CC  = \
	gauss2.c  gauss3_hcn.c	gauss7_n2h+.c  gaussn.c  myline.c  taugauss.c rotcurm.c 

FIT_CXX = \

FIT_F77 = \


SRCFILES= Makefile \
	$(FIT_CC) $(FIT_F77) $(FIT_CXX)

O_FILES= $(FIT_CC:c=o) $(FIT_CXX:cc=o) $(FIT_F77:f=o)

OBJFILES= $(O_FILES:o=so)

BINFILES = 

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`

clean:
	rm -f *.o *.so *.a core $(BINFILES) $(TESTFILES) ctextest.*

# all:	

obj:	$(OBJFILES)

#	update the library
update: .update_lib

.update_lib: $(OBJFILES)
	mv $? $(NEMOOBJ)/fit
	@touch .update_lib

#	install the library: can be (re)done from scratch by removing .o files
install:  .install_lib

.install_lib: $(OBJFILES)
	cp $? $(NEMOOBJ)/fit

.install_src: $(SRCFILES)
	cp $? $(NEMO)/src/kernel/tab/fit
	@touch .install_src

# Targets used by Nemo to update libraries and bins.

all:

nemo_lib: $(LIBFILE)


nemo_sublibs:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_lib); done

nemo_bin: nemo_subbins $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	rm -f *.o

bins:

nemo_subbins:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_bin); done

nemo_src:
	-@for i in $(SUBDIRS); do \
	(cd $$i; $(MAKE) nemo_src); done
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

#

tests: 

