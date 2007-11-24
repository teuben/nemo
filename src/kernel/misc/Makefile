#	Makefile for NEMO/src/kernel/misc
include $(NEMOLIB)/makedefs

#	NUNEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

LPP = $(NEMOLIB)/libnemo++.a
.PRECIOUS: $(LPP)

# Files etc. in this directory

MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = axis.h hash.h vectmath.h cgs.h mks.h layout.h
SRCFILES= axis.c besselfunc.c erf.c fie.c \
	  frandom.c grid.c \
	  hash.c herinp.c layout.c linreg.c log2.c \
	  lsq.c matinv.c nemofie.c imsl.c \
	  match.c mdarray.c median.c minmax.c moment.c \
	  nemoinp.c nemomain.c newextn.c pick.c pow.c scanopt.c \
	  setfblank.c spline.c sqr.c timers.c vectmath.c within.c \
	  xrand.c xrandom.c \
	  sort.c sortptr.c unwrap.c

OBJFILES= axis.o besselfunc.o erf.o fie.o \
	  frandom.o grid.o \
	  hash.o herinp.o layout.o linreg.o log2.o \
	  lsq.o matinv.o nemofie.o imsl.o \
	  match.o mdarray.o median.o minmax.o moment.o \
	  nemoinp.o nemomain.o newextn.o pick.o pow.o scanopt.o \
	  setfblank.o spline.o sqr.o timers.o vectmath.o within.o \
	  xrand.o xrandom.o \
	  sort.o sortptr.o unwrap.o

LOBJFILES= $L(axis.o) $L(besselfunc.o) $L(erf.o) $L(fie.o) $L(layout.o) \
	  $L(frandom.o) $L(grid.o) \
	  $L(hash.o) $L(herinp.o) $L(linreg.o) $L(log2.o) \
	  $L(lsq.o) $L(matinv.o) $L(nemofie.o) $L(imsl.o) \
	  $L(match.o) $L(mdarray) $L(median.o) $L(minmax.o) $L(moment.o) \
	  $L(nemoinp.o) $L(nemomain.o) $L(newextn.o) $L(pick.o) $L(pow.o) $L(scanopt.o) \
	  $L(setfblank.o) $L(spline.o) $L(sqr.o) $L(timers.o) $L(vectmath.o) $L(within.o) \
	  $L(xrand.o) $L(xrandom.o) \
	  $L(sort.o) $L(sortptr.o) $L(unwrap.o)

BINFILES = nemoinp nemoshow layout xrandom scanopt 

TESTFILES = vecttest axistest splinetest withintest \
	matchtest linreg momenttest gridtest unwraptest frandomtest \
	mdarraytest timerstest

#	update the library: direct comparison with modules inside L
help:
	@echo Nemo makefile in `pwd`
	@echo bin: BINFILES=$(BINFILES)
	@echo tests: TESTFILES=$(TESTFILES)
	@echo install: OBJFILES=$(OBJFILES)

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
install:   .install_lib   .install_lib++

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_lib++: nemomain++.o
	ar ruv $(LPP) $?
	$(RANLIB) $(LPP)

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)
	@touch .install_inc

.install_bin: $(BINFILES)
	mv $? $(NEMOBIN)

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

vectmath_d.o: vectmath.h vectmath.c
	$(CC) $(CFLAGS) -USINGLEPREC -c vectmath.c
	mv vectmath.o vectmath_d.o

vectmath_f.o: vectmath.h vectmath.c
	$(CC) $(CFLAGS) -DSINGLEPREC -c vectmath.c
	mv vectmath.o vectmath_f.o

# Nostandard targets

nemoinp: nemoinp.c
	$(CC) $(CFLAGS) -o nemoinp nemoinp.c -DTOOLBOX  $(NEMO_LIBS) -lm

herinp: 
	$(CC) $(CFLAGS) -o herinp herinp.c -DTESTBED -DNEED_DCD $(NEMO_LIBS) -lm

layout: layout.c
	$(CC) $(CFLAGS) -o layout layout.c -DTOOLBOX $(NEMO_LIBS) $(YAPPLIB) -lm

obj: $(OBJFILES)
bin: $(BINFILES)
tests: $(TESTFILES) 

vecttest: vecttest.c vectmath.h vectmath.c
	$(CC) $(CFLAGS) -o vecttest vecttest.c  $(NEMO_LIBS)

xrandom: xrandom.c 
	$(CC) $(CFLAGS) -o xrandom -DTOOLBOX xrandom.c $(NEMO_LIBS)

scanopt: scanopt.c
	$(CC) $(CFLAGS) -o scanopt -DTOOLBOX scanopt.c $(NEMO_LIBS)

xrandomtest: xrandom.c 
	$(CC) $(CFLAGS) -o xrandomtest -DTESTBED xrandom.c $(NEMO_LIBS)

frandomtest: frandom.c 
	$(CC) $(CFLAGS) -o frandomtest -DTESTBED frandom.c $(NEMO_LIBS)

hashtest: hash.c 
	$(CC) $(CFLAGS) -o hashtest -DTESTBED hash.c $(NEMO_LIBS)

mdarraytest: mdarray.c 
	$(CC) $(CFLAGS) -o mdarraytest -DTESTBED mdarray.c $(NEMO_LIBS)

axistest: axis.c
	$(CC) $(CFLAGS) -o axistest -DTESTBED axis.c $(NEMO_LIBS)

axistest_ps: axis.c
	$(CC) $(CFLAGS) -o axistest_ps -DTESTBED axis.c $(NEMO_LIBS)

splinetest: spline.c
	$(CC) $(CFLAGS) -o splinetest -DTESTBED spline.c $(NEMO_LIBS)

withintest: within.c
	$(CC) $(CFLAGS) -o withintest -DTESTBED within.c $(NEMO_LIBS)

besseltest: besselfunc.c
	$(CC) $(CFLAGS) -o besseltest -DTESTBED besselfunc.c $(NEMO_LIBS)

matchtest: match.c
	$(CC) $(CFLAGS) -o matchtest -DTESTBED match.c $(NEMO_LIBS)

linreg: linreg.c
	$(CC) $(CFLAGS) -o linreg -DTESTBED linreg.c $(NEMO_LIBS)

momenttest: moment.c
	$(CC) $(CFLAGS) -o momenttest -DTESTBED moment.c $(NEMO_LIBS)

gridtest: grid.c
	$(CC) $(CFLAGS) -o gridtest -DTESTBED grid.c $(NEMO_LIBS)

unwraptest: unwrap.c
	$(CC) $(CFLAGS) -o unwraptest -DTESTBED unwrap.c $(NEMO_LIBS)

timerstest: timers.c
	$(CC) $(CFLAGS) -o timerstest -DTESTBED timers.c $(NEMO_LIBS)

powtest: pow.c
	$(CC) $(CFLAGS) -o powtest -DTESTBED pow.c $(NEMO_LIBS)


# special targets

bad_fortran:
	make install OBJFILES=nemomain.o CFLAGS="-g -DMAIN"
