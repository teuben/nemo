#	Makefile for NEMO/src/orbit/potential/data
include $(NEMOLIB)/makedefs

#	Compile flags  $(NFLAGS)
DIM=THREEDIM
NFLAGS=-D$(DIM)
LOCAL_INC = $(NFLAGS) -Iinc -fpic

#	remote dir where object files should live
R = $(NEMOOBJ)/potential

#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# testing a single potential
TEST = hubble.c

# Files etc. in this directory

MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 

POT_CC = \
	bar83.c brandt.c bulge1.c ccd.c cp80.c dehnen.c dublinz.c expdisk.c flatz.c halo.c \
	grow_plum.c grow_plum2.c hdgrow1.c \
	harmonic.c hernquist.c hh64.c hom.c hubble.c kim11.f kuzmindisk.c \
	isochrone.c jaffe.c log.c log2.c mestel.c miyamoto.c nfw.c nfw2.c null.c \
	op73.c plummer.c plummer2.c persic.c rh84.c rotcur0.c rotcur1.c rotcure.c rotcurm.c rotcur.c \
	sh76.c teusan85.c triax.c \
	turner92.c twofixed.c plummer4.c vertdisk.c tidaldisk.c polynomial.c wada94.c \
        zero.c point.c

POT_F77 = \
	athan92.f harmonicf.f pfenniger84.f piner94.f 

SRCFILES= Makefile \
	$(POT_CC) $(POT_F77) $(POT_CXX)

O_FILES= $(POT_CC:c=o) $(POT_CXX:cc=o) $(POT_F77:f=o)

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
	mv $? $(NEMOOBJ)/potential
	@touch .update_lib

#	install the library: can be (re)done from scratch by removing .o files
install:  .install_lib

.install_lib: $(OBJFILES)
	cp $? $(NEMOOBJ)/potential

.install_src: $(SRCFILES)
	cp $? $(NEMO)/src/orbit/potential/data
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

potctex.inc: $(POT_CC)
	$(NEMO)/src/scripts/ctex -latex -section -nohead -nofoot $(POT_CC) > potctex.inc
	-diff potctex.inc $(NEMO)/text/manuals

ctex:
	$(NEMO)/src/scripts/ctex -latex -section $(TEST) > ctextest.tex
	latex ctextest.tex
	dvips ctextest -o
	@echo Test for TEST=$(TEST) done... results in ctextest
