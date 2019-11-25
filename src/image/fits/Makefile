#	Makefile for NEMO/src/image/fits
include $(NEMOLIB)/makedefs

LOCAL_LIB = $(CFITSIO_LIB)

#PRECFLAG = -DSINGLEPREC -fsingle
PRECFLAG = 

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = fitsio.h fits.h
SRCFILES = fitsio.c fits.c \
	   ccdfits.c fitsccd.c fitssplit.c scanfits.c
OBJFILES=  fitsio_nemo.o fits.o
LOBJFILES= $L(fitsio.o) $L(fits.o)
BINFILES = ccdfits fitsccd fitssplit scanfits fitstab fitshead \
	   fitsglue fits8to16 tabfits
# fitsgrid <-- is broken ?
TESTFILES=  format

ifneq ($(strip $(HDF_LIB)),)
  BINFILES += sdsfits 
endif

help:
	@echo NEMO V2.2 NEMO/src/image/fits
	@echo bin: BINFILES=$(BINFILES)
	@echo test: TESTFILES=$(TESTFILES)

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



# special

fitsgids: fitsgids.c look_nemo.c
	$(CC) $(CFLAGS) -o fitsgids fitsgids.c look_nemo.c \
		 $(NEMO_LIBS) -lgdi

sdsfits: sdsfits.c 
	$(CC) $(CFLAGS) -o sdsfits sdsfits.c $(HDF_LIB) $(NEMO_LIBS)  $(CFITSIO_LIB) $(EL)

mkfits: mkfits.c 
	$(CC) $(CFLAGS) -o mkfits mkfits.c $(NEMO_LIBS) $(CFITSIO_LIB) $(EL)

sdinfo: sdinfo.c
	$(CC) $(CFLAGS) -o sdinfo sdinfo.c $(NEMO_LIBS) $(CFITSIO_LIB) $(EL)

ccdfitsio: ccdfitsio.c 
	$(CC) $(CFLAGS) -o ccdfitsio ccdfitsio.c $(NEMO_LIBS) $(CFITSIO_LIB) $(EL)


fitstiff: fitstiff.c 
	$(CC) $(CFLAGS) -o fitstiff fitstiff.c \
		 $(NEMO_LIBS) $(EL) -ltiff


ccd_cfitsio.o:		ccd_cfitsio.c
	$(CC) $(OPT) $(CFITSIOINC) -c $*.c


fitsio_nemo.o:	fitsio.c fitsio_nemo.c
