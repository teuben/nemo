#	The FORTRAN interface routines
include $(NEMOLIB)/makedefs

L = $(NEMOLIB)/libnemo.a

SRCS = Makefile README ftoc.c getparam_ftoc.c \
	test.f \
	getparam.inc
INCS = getparam.inc
BINS = ftoc test unfio
OBJFILES = getparam_ftoc.o unfio.o
TESTFILES = testunfio funcc funcf

SRCDIR = $(NEMO)/src/kernel/fortran

help:
	@echo "Testing still out here - these routines are not "
	@echo "Added to the library by default - use 'make install'"
	@echo tests: TESTFILES=$(TESTFILES)
	@echo bins: BINS=$(BINS)


install: .install_lib


.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_src: $(SRCS)
	cp $? $(SRCDIR)
	@touch .install_src

.install_inc: $(INCS)
	cp $? $(NEMOINC)
	@touch .install_inc

clean:
	rm -f *~ *.o *.a core $(BINS) $(TESTFILES) test_main.c mainfuncf_tmp.*

bins:	$(BINS)
tests:  $(TESTFILES)
	

test:	test.f getparam_ftoc.o ftoc
	ftoc test.f > test_main.c
	$(CC) $(CFLAGS) -c test_main.c
	$(FC) $(FFLAGS) -c test.f
	$(FC) $(FFLAGS) -o test test.o getparam_ftoc.o test_main.o -L$(NEMOLIB) -lnemo
	nm test.o | grep -i nemomain
	./test test.f

ftoc:	ftoc.c
	$(CC) $(CFLAGS) -o ftoc ftoc.c \
	  -L$(NEMOLIB) -lnemo -lm

# getparam_ftoc.o: getparam_ftoc.c



unfio:	unfio.c
	$(CC) -DTOOLBOX -o unfio unfio.c -lnemo -lm

testunfio: testunfio.f
	$(FC) $(FFLAGS) -o testunfio testunfio.f


#	C main calls a Fortran routine ; use FLIBS to link
funcc:
	$(CC) $(CFLAGS) -c funcc.c
	$(FC) $(FFLAGS) -c funcf.f
	$(CC) $(CFLAGS) -c mainfuncc.c
	$(CC) $(LDFLAGS) -o funcc mainfuncc.o funcc.o funcf.o $(NEMO_LIBS) $(FLIBS)

funcf:
	$(CC) $(CFLAGS) -c funcc.c
	$(FC) $(FFLAGS) -c funcf.f
	$(FC) $(FFLAGS) -c mainfuncf.f
	ftoc mainfuncf.f > mainfuncf_tmp.c
	$(CC) $(CCFLAGS) -c mainfuncf_tmp.c
	$(CC) $(LDFLAGS) -o funcf mainfuncf_tmp.o mainfuncf.o funcc.o funcf.o $(NEMO_LIBS) $(FLIBS)


