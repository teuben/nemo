#	Makefile for NEMO/src/image/fits
include $(NEMOLIB)/makedefs

LOCAL_LIB = -lnetcdf

#PRECFLAG = -DSINGLEPREC -fsingle
PRECFLAG = 

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = lmtinfo.c
OBJFILES=  
LOBJFILES= 
BINFILES = lmtinfo
TESTFILES=  

help:
	@echo NEMO V2.2 NEMO/src/image/cdf
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

lmtinfo : lmtinfo.c 
	$(CC) $(CFLAGS) -o lmtinfo lmtinfo.c \
		$(NEMO_LIBS) $(LOCAL_LIB) $(EL)


