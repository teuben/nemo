#	Makefile for NEMO/src/kernel/numrec (subset of Numerical Recipes)
include $(NEMOLIB)/makedefs
#	some routines exist in single as well as double precision


#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory

MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = nrutil.h nr.h
SRCFILES= eigsrt.c   jacobi.c nrutil.c \
	  eigsrt_d.c jacobi_d.c ran3.c beta.c gamma.c pearsn.c \
	  polint.c trapzd.c qromb.c \
	  gaussj.c mrqcof.c covsrt.c mrqmin.c lfit.c \
	nr_nllsqfit.c



OBJFILES= eigsrt.o jacobi.o nrutil.o \
	  eigsrt_d.o jacobi_d.o ran3.o beta.o gamma.o pearsn.o \
	  probks.o kstwo.o \
	  polint.o trapzd.o qromb.o \
	  gaussj.o mrqcof.o covsrt.o mrqmin.o lfit.o \
	nr_nllsqfit.o

LOBJFILES= $L(eigsrt.o)   $L(jacobi.o) $L(nrutil.o) \
	   $L(eigsrt_d.o) $L(jacobi_d.o) $L(ran3.o) \
	   $L(beta.o) $L(gamma.o) $L(pearsn.o) \
	  $L(polint.o) $L(trapzd.o) $L(qromb.o)

BINFILES = 
TESTFILES =

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`


clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

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
install:  .install_lib 

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)

bins:

# Targets used by Nemo to update libraries and bins.

all:

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

vectmath_d.o: vectmath.h vectmath.c
	$(CC) $(CFLAGS) -USINGLEPREC -c vectmath.c
	mv vectmath.o vectmath_d.o

vectmath_f.o: vectmath.h vectmath.c
	$(CC) $(CFLAGS) -DSINGLEPREC -c vectmath.c
	mv vectmath.o vectmath_f.o

# Nostandard targets

nemoinp: nemoinp.c
	$(CC) $(CFLAGS) -o nemoinp nemoinp.c -DTOOLBOX  $(NEMO_LIBS)

tests: $(TESTFILES) 

vecttest: vecttest.c vectmath.h vectmath.c
	$(CC) $(CFLAGS) -o vecttest vecttest.c  $(NEMO_LIBS)

xrandomtest: xrandom.c 
	$(CC) $(CFLAGS) -o xrandomtest -DTESTBED xrandom.c  $(NEMO_LIBS)

frandomtest: frandom.c 
	$(CC) $(CFLAGS) -o frandomtest -DTESTBED frandom.c  $(NEMO_LIBS)

axistest: axis.c
	$(CC) $(CFLAGS) -o axistest -DTESTBED axis.c $(NEMO_LIBS)

axistest_ps: axis.c
	$(CC) $(CFLAGS) -o axistest_ps -DTESTBED axis.c  $(NEMO_LIBS)

splinetest: spline.c
	$(CC) $(CFLAGS) -o splinetest -DTESTBED spline.c $(NEMO_LIBS)

withintest: within.c
	$(CC) $(CFLAGS) -o withintest -DTESTBED within.c  $(NEMO_LIBS)

scanopttest: scanopt.c
	$(CC) $(CFLAGS) -o scanopttest -DTESTBED scanopt.c  $(NEMO_LIBS)

besseltest: besselfunc.c
	$(CC) $(CFLAGS) -o besseltest -DTESTBED besselfunc.c  $(NEMO_LIBS)

