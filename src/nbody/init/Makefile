#	Makefile for NEMO/src/nbody/init
include $(NEMOLIB)/makedefs
# LOCAL_LIB = $(FORLIBS)

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = 
OBJFILES=  
LOBJFILES= 
BINFILES = mkplummer mkdisk mkconfig mkexpdisk \
	   mkbaredisk mkommod mkpolytrope mktestdisk \
	   plummer mkexphot snapenter mkspiral mkhom mkhomsph \
	   mkop73 mkcube mksphere mkflowdisk magalie
TESTFILES= 

help:
	@echo NEMO $$NEMO/src/nbody/init

clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

all:	$(BINFILES) $(TESTFILES)

bin:	$(BINFILES)

#	update or install a library
update:  .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

install:

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
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

#	Special rules

mkplummer: mkplummer.c
	$(CC) $(CFLAGS) -DTOOLBOX -o mkplummer mkplummer.c \
			$(NEMO_LIBS) $(EL) -lm

mkkd95:	kd95

magalie: bkpg2001

kd95:
	(cd $(NEMO)/usr/kuijken/GalactICS-exp/src; make -f Makefile.nemo clean all install)

bkpg2001:
	(cd $(NEMO)/usr/boily/magalie ; $(MAKE) install)
