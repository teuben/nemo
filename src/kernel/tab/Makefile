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
SRCFILES = table.c gettab.c tabselect.c funtab.c getaline.c
OBJFILES=  table.o gettab.o tabselect.o funtab.o getaline.o
LOBJFILES= $L(table.o) $L(gettab.o) $L(tabselect.o) $L(funtab.o) $L(getaline.o)
BINFILES = tabhist tablst tabplot tablsqfit tabmath gettab funtab \
	   tabcomment tabspline
TESTFILES= 

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
		$(L) $(YAPPLIB) $(EL) -lm

tabview: tabview.c
	$(CC) $(CFLAGS) -o tabview tabview.c \
		$(L) $(YAPPLIB) $(EL) -lm

tabhist: tabhist.c
	$(CC) $(CFLAGS) -DYAPP -o tabhist tabhist.c \
		$(L) $(YAPPLIB) $(EL) -lm

tablsqfit: tablsqfit.c
	$(CC) $(CFLAGS) -o tablsqfit tablsqfit.c \
		$(L) $(EL) -lm

tabspline: tabspline.c
	$(CC) $(CFLAGS) -o tabspline tabspline.c \
		$(L) $(EL) -lm

gettab: gettab.c
	$(CC) -DTOOLBOX $(CFLAGS) -o gettab gettab.c \
		$(L) $(EL) -lm

funtab: funtab.c
	$(CC) -DTOOLBOX $(CFLAGS) -o funtab funtab.c \
		$(L) $(EL) -lm



