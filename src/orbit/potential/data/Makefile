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
SRCFILES= Makefile \
	athan92.f \
	bar83.c bulge1.c ccd.c cp80.c dehnen.c dublinz.c expdisk.c flatz.c halo.c \
	Combined.cc Dehnen.cc Hernquist.cc NFW.cc GalPot.cc MiyamotoNagai.cc \
	Plummer.cc \
	grow_plum.c grow_plum2.c \
	harmonic.c harmonicf.f hernquist.c hom.c hubble.c kuzmindisk.c \
	isochrone.c jaffe.c log.c mestel.c miyamoto.c nfw.c null.c \
	op73.c pfenniger84.f piner94.f plummer.c plummer2.c rh84.c rotcur0.c rotcur.c teusan85.c triax.c \
	twofixed.c plummer4.c vertdisk.c tidaldisk.c polynomial.c wada94.c
OBJFILES=\
	athan92.so \
	bar83.so bulge1.so ccd.so cp80.so dehnen.so dublinz.so expdisk.so flatz.so halo.so \
	Combined.so \
	Dehnen.so Hernquist.so NFW.so GalPot.so MiyamotoNagai.so Plummer.so \
	grow_plum.so grow_plum2.so \
	harmonic.so harmonicf.so hernquist.so hom.so hubble.so kuzmindisk.so \
	isochrone.so jaffe.so log.so mestel.so miyamoto.so nfw.so null.so \
	op73.so pfenniger84.o piner94.so plummer.so plummer2.so rh84.so rotcur0.so rotcur.so teusan85.so triax.so \
	twofixed.so plummer4.so vertdisk.so tidaldisk.so polynomial.so wada94.so

BINFILES = 

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`

clean:
	rm -f *.o *.so *.a core $(BINFILES) $(TESTFILES) ctextest.*

all:	

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


ctex:
	$(NEMO)/src/scripts/ctex -latex -section $(TEST) > ctextest.tex
	latex ctextest.tex
	dvips ctextest -o
	@echo Test for TEST=$(TEST) done... results in ctextest
