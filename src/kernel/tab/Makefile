#	Makefile for NEMO/src/kernel/tab
include $(NEMOLIB)/makedefs

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = 
SRCFILES = table.c gettab.c tabselect.c funtab.c getaline.c pyplot.c
OBJFILES=  table.o gettab.o tabselect.o funtab.o getaline.o pyplot.o
LOBJFILES= $L(table.o) $L(gettab.o) $L(tabselect.o) $L(funtab.o) $L(getaline.o) $L(pyplot.o)
BINFILES = tabhist tablst tabplot tablsqfit tabmath gettab funtab meanmed \
	   tabcomment tabspline tab2xml tabnllsqfit tabdate tabfilter tabtrend \
	   tabstat
TESTFILES= getaline

help:
	@echo NEMO/src/kernel/tab

clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

#	update or install a library
update: .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

$(L):	$(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?

install:  .install_lib 

.install_lib: $(OBJFILES) 
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)
	@touch .install_inc

bins:	$(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

# Targets used by Nemo to update libraries and bins.

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

bin: $(BINFILES)
tests: $(TESTFILES)

# special

tabplot: tabplot.c
	$(CC) $(CFLAGS) -o tabplot tabplot.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

tabplot_command: tabplot_command.c
	$(CC) $(CFLAGS) -o tabplot_command tabplot_command.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

tabview: tabview.c
	$(CC) $(CFLAGS) -o tabview tabview.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

tabzoom: tabzoom.c
	$(CC) $(CFLAGS) -o tabzoom tabzoom.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

tabhist: tabhist.c
	$(CC) $(CFLAGS) -DYAPP -o tabhist tabhist.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lm

tabhist_flogger: tabhist.c
	$(CC) $(CFLAGS) -DYAPP -DFLOGGER -o tabhist_flogger tabhist.c \
		$(NEMO_LIBS) $(YAPPLIB) $(EL) -lsort -lm

tablsqfit: tablsqfit.c
	$(CC) -UHAVE_GSL $(CFLAGS) -o tablsqfit tablsqfit.c \
		$(NEMO_LIBS) $(EL) -lm

MY_GSL = -lgsl -lgslcblas 

tablsqfit_gsl: tablsqfit.c
	$(CC) -DHAVE_GSL $(CFLAGS) -o tablsqfit_gsl tablsqfit.c \
		$(NEMO_LIBS) $(GSL_LIBS) $(EL) -lm

tabnllsqfit: tabnllsqfit.c
	$(CC) $(CFLAGS) -o tabnllsqfit tabnllsqfit.c \
		$(NEMO_LIBS) $(EL) -lm

tabspline: tabspline.c
	$(CC) $(CFLAGS) -o tabspline tabspline.c \
		$(NEMO_LIBS) $(EL) -lm

tab2xml: tab2xml.c
	$(CC) $(CFLAGS) -o  tab2xml tab2xml.c \
		$(NEMO_LIBS) $(EL) -lm

tabpoly: tabpoly.c
	$(CC) $(CFLAGS) -o tabpoly tabpoly.c \
		$(NEMO_LIBS) $(EL) -lm

gettab: gettab.c
	$(CC) -DTOOLBOX $(CFLAGS) -o gettab gettab.c \
		$(NEMO_LIBS) $(EL) -lm

funtab: funtab.c
	$(CC) -DTOOLBOX $(CFLAGS) -o funtab funtab.c \
		$(NEMO_LIBS) $(EL) -lm

tabfft: tabfft.c
	$(CC) $(CFLAGS) -o tabfft tabfft.c \
		$(NEMO_LIBS) $(EL) -lm


getaline:	getaline.c
	$(CC) -DTESTBED $(CFLAGS) -o getaline getaline.c \
		$(NEMO_LIBS) $(EL) -lm

