#	Makefile for NEMO/src/image/io
include $(NEMOLIB)/makedefs

#PRECFLAG = -DSINGLEPREC -fsingle
PRECFLAG = 

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = image.h matdef.h xyio.h
SRCFILES = image.c ccddump.c xyio.c 
OBJFILES=  image.o xyio.o wcsio.o
LOBJFILES= $L(image.o) $L(xyio.o) $L(wcsio.o)
BINFILES = ccddump ccdprint ccdslice sigccd
TESTFILES= imagetest xyiotest

help:
	@echo NEMO V2.2 NEMO/src/image/io
	@echo OBJFILES: $(OBJFILES)
	@echo BINFILES: $(BINFILES)
clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)
bin: $(BINFILES)
test: $(TESTFILES)

#	update or install a library
update:  .update_inc .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

.install_bin: $(BINFILES)
	mv $? $(NEMOBIN)

$(L):	$(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

install:   .install_lib 

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

#	Special rules

tests: $(TESTFILES)

imagetest:  image.c
	$(CC) $(CFLAGS) -DTESTBED -o imagetest image.c $(L) -lm

xyiotest:  xyio.c
	$(CC) $(CFLAGS) -DTESTBED -o xyiotest xyio.c $(L) -lm

# special

