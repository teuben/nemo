#	Makefile for NEMO/src/kernel/cores
include $(NEMOLIB)/makedefs

#	Compile flags
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
INCFILES = filefn.h options.h stdinc.h strlib.h 
SRCFILES= allocate.c  bswap.c burststring.c common.c cputime.c date_id.c \
	  error.c file_size.c filefn.c strlib.c progress.c
OBJFILES= allocate.o bswap.o burststring.o common.o cputime.o date_id.o \
	  error.o file_size.o filefn.o strlib.o progress.o
##  exception.o		
# exception doesn't work on solaris
LOBJFILES= $L(allocate.o)  $L(burststring.o) $L(cputime.o) $L(date_id.o) \
	   $L(error.o) $L(file_size.o) $L(filefn.o) $L(strlib.o) $L(bswap.o) \
	   $L(progress.o)
BINFILES = bswap
TESTFILES= errortest burstrtest filefntest sizetest cputest cputool \
	   progresstest allocatetest

help:
	@echo NEMO/src/kernel/cores
	@echo bins: BINFILES=$(BINFILES)
	@echo tests: TESTFILES=$(TESTFILES)

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

# Targets used to assemble object library.
bin: $(BINFILES)
obj: $(OBJFILES)
tests: $(TESTFILES)

allocatetest: allocate.c
	$(CC) $(CFLAGS) -o allocatetest -DTESTBED allocate.c $(NEMO_LIBS)

errortest: error.c
	$(CC) $(CFLAGS) -o errortest -DTESTBED error.c $(NEMO_LIBS)

progresstest: progress.c
	$(CC) $(CFLAGS) -o progresstest -DTESTBED progress.c $(NEMO_LIBS)

burstrtest: burststring.c
	$(CC) $(CFLAGS) -o burstrtest -DTESTBED burststring.c $(NEMO_LIBS)

filefntest: filefn.c
	$(CC) $(CFLAGS) -o filefntest -DTESTBED filefn.c $(NEMO_LIBS)

swaptest: bswap.c
	$(CC) $(CFLAGS) -o swaptest -DTESTBED bswap.c $(NEMO_LIBS)

bswap:	bswap.c
	$(CC) $(CFLAGS) -o bswap -DTOOLBOX bswap.c $(NEMO_LIBS)

sizetest: file_size.c
	$(CC) $(CFLAGS) -o sizetest -DTESTBED file_size.c $(NEMO_LIBS)

# cputest is old nemo testbed
cputest: cputime.c
	$(CC) $(CFLAGS) -o cputest -DTESTBED cputime.c $(NEMO_LIBS)

# cputool is starlab toolbox
cputool: cputime.c
	$(CC) $(CFLAGS) -o cputool -DTOOLBOX cputime.c $(NEMO_LIBS)

