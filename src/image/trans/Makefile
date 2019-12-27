#	Makefile for NEMO/src/image/trans
include $(NEMOLIB)/makedefs

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
BINFILES = ccdsmooth ccdmath ccdflip ccdfill ccdsharp ccdsharp3 \
	ccdclip ccdintpol ccdmedian ccdpot ccdgen
TESTFILES= 

help:
	@echo NEMO V2.2 NEMO/src/image/trans

clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

bin: $(BINFILES)

#	update or install a library
update:  .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
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

.install_bin: $(BINFILES)
	mv $? $(NEMOBIN)

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


ccdplot:  ccdplot.c
	$(CC) $(CFLAGS) -o ccdplot ccdplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

ccdgoat:  ccdgoat.c
	$(CC) $(CFLAGS) -o ccdgoat ccdgoat.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm


bench1:
	$(TIME) ccdreorder . . dims=100,100,100  iter=512 order=0,1,2 help=cm
	$(TIME) ccdreorder . . dims=100,100,100  iter=512 order=2,1,0 help=cm
	$(TIME) ccdreorder . . dims=200,200,200  iter=64  order=0,1,2 help=cm
	$(TIME) ccdreorder . . dims=200,200,200  iter=64  order=2,1,0 help=cm
	$(TIME) ccdreorder . . dims=400,400,400  iter=8   order=0,1,2 help=cm
	$(TIME) ccdreorder . . dims=400,400,400  iter=8   order=2,1,0 help=cm
	$(TIME) ccdreorder . . dims=800,800,800  iter=1   order=0,1,2 help=cm
	$(TIME) ccdreorder . . dims=800,800,800  iter=1   order=2,1,0 help=cm
