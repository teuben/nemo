#	Makefile for NEMO/src/nbody/reduc
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
BINFILES = snapplot snapplot3 snapdiagplot snapmstat snapprint \
    snapcmp snapcmphist snapcmpplot snapsample snapstat \
    snapmnmx radprof snapkinplot \
    snapfour snapkinem snapmradii snapopt snaprstat snaptrak  \
    snapvratio snapshell snapplotv snapfit snapkmean
TESTFILES= 

help:
	@echo NEMO V2.2 NEMO/src/nbody/reduc

clean:
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

all:	$(BINFILES) $(TESTFILES)
bin: $(BINFILES)
test: $(TESTFILES)
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

trakplot:
	(cd $(NEMOBIN); ln -s snapplot trakplot)

snapplot:  snapplot.c
	$(CC) $(CFLAGS) -o snapplot snapplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapplot3:  snapplot3.c
	$(CC) $(CFLAGS) -o snapplot3 snapplot3.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapplot_cg:  snapplot.c
	$(CC) -DCOLOR $(CFLAGS) -o snapplot_cg snapplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

aasplot:  aasplot.c
	$(CC) -DCOLOR $(CFLAGS) -o aasplot aasplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapplotv:  snapplotv.c
	$(CC) $(CFLAGS) -o snapplotv snapplotv.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapdiagplot:  snapdiagplot.c
	$(CC) $(CFLAGS) -o snapdiagplot snapdiagplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

radprof:  radprof.c
	$(CC) $(CFLAGS) -o radprof radprof.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapkinplot:  snapkinplot.c
	$(CC) $(CFLAGS) -o snapkinplot snapkinplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapinert:  snapinert.c eigen.c
	$(CC) $(CFLAGS) -I$(NUMRECPATH) -o snapinert snapinert.c eigen.c \
		$(NEMO_LIBS) $(NUMRECLIB) $(EL) -lm

snapcmp:  snapcmp.c
	$(CC) $(CFLAGS) -DSTATISTICS -o snapcmp snapcmp.c \
		$(NEMO_LIBS) $(EL) -lm

snapcmphist:  snapcmp.c
	$(CC) $(CFLAGS) -DHISTOGRAM -o snapcmphist snapcmp.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

snapcmpplot:  snapcmp.c
	$(CC) $(CFLAGS) -DSCATTERPLOT -o snapcmpplot snapcmp.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

