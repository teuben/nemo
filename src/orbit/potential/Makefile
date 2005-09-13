#	Makefile for NEMO/src/orbit/potential
include $(NEMOLIB)/makedefs
LOCAL_LIB = $(FORLIBS)

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# local for potentials (won't work well on SGI if .o is actually shared object)
EL = 
POTNAME = zero
POTOBJ = $(NEMOOBJ)/potential/$(POTNAME).so
#POTOBJ =

# Files etc. in this directory

MAN1FILES = potlist.1 rotcurves.1
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES= potential.c potential0.c potlist.c rotcurves.c acceleration.cc \
          $(INCFILES) Makefile
OBJFILES= potential.o fmath.o acceleration.o
LOBJFILES= $L(potential.o) $L(fmath.o)
BINFILES = potlist rotcurves potccd potq potrot potsf
TESTFILES = potfortest

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`
	@echo Example:
	@echo '  make potlist POTOBJ="data/athan92.o $$NEMOLIB/potential0.o"'
	@echo "NOTE: for some architectures (e.g. dec alpha)"
	@echo "linking with a shared object POTOBJ"
	@echo "	POTOBJ=$(POTOBJ)"
	@echo "should not be done: all programs need 'make ... POTOBJ= '"
	@echo "set to blank"


clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

all:	$(BINFILES) $(TESTFILES)
bin:	$(BINFILES)

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
install:   .install_lib 

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_inc:

.install_src: $(SRCFILES)
	cp $? $(NEMO)/src/orbit/potential
	@touch .install_src


###	$(CC) $(CFLAGS) -o $* $*.c $(POTOBJ) $(L) $(EL) $(FORLIBS) -lm

#
# Targets used to export things to Nemo.

#
# Targets used by Nemo to update libraries and bins.

all:

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

# Targets used to assemble object library.

rotcurves: rotcurves.c
	$(CC) $(CFLAGS) -o rotcurves rotcurves.c $(NEMO_LIBS) $(EL) $(YAPPLIB) $(FORLIBS) -lm

rotcurves_ps: rotcurves.c
	$(CC) $(CFLAGS) -o rotcurves_ps rotcurves.c $(NEMO_LIBS) $(EL) $(YAPPLIB_PS) $(FORLIBS) -lm
#

tests: 
        

potential.src:
	ixfile potential.{3,5,h} > potential.src
	cat potential.c >> potential.src

# use POTNAME= from data

potfortest: potentialf.o
	@echo POTNAME=$(POTNAME)
	$(FC) $(FFLAGS) -o potfortest potfortest.f potentialf.o $(POTOBJ)

# 
#  C++ testing; these are clones of e.g. potlist, but for double and float
#
try1: try1.cc
	$(CXX) $(CFLAGS) -o try1 try1.cc $(NEMO_LIBSPP) $(EL) $(FORLIBS) -lm

try2: try2.cc
	$(CXX) $(CFLAGS) -o try2 try2.cc $(NEMO_LIBSPP) $(EL) $(FORLIBS) -lm


potlist++: potlist.c
	$(CXX) $(CFLAGS) -o potlist++ potlist.c $(NEMO_LIBSPP) $(EL) $(FORLIBS) -lm
