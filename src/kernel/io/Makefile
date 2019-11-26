#	Makefile for NEMO/src/kernel/io
include $(NEMOLIB)/makedefs

#	sub-directories to be visited, in the order given here
#DIRS = 

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files etc. in this directory
MAN1FILES = 
MAN3FILES = 
MAN5FILES = 
INCFILES = story.h
SRCFILES = dprintf.c command.c convert.c cvsid.c defv.c endian.c extstring.c \
	   filesecret.[ch] getparam.[ch] history.[ch] memio.c outdefv.c \
	   story.[ch] stropen.c mstropen.c usage.c \
	   ieeehalfprecision.c \
	   filestruct.h Makefile
OBJFILES=  dprintf.o command.o convert.o cvsid.o defv.o endian.o extstring.o \
	   filesecret.o getparam.o history.o memio.o outdefv.o \
	   ieeehalfprecision.o \
	   stropen.o mstropen.o usage.o 
LOBJFILES= $L(dprintf.o) $L(command.o) $L(convert.o) $L(cvsid.o) $L(defv.o) $L(endian.o) $L(extstring.o) \
           $L(filesecret.o) $L(getparam.o) $L(history.o) $L(memio.o) $L(outdefv.o) \
	   $L(ieeehalfprecision.o) $L(stropen.o) $L(mstropen(.o) $L(usage.o)
BINFILES = csf tsf rsf qsf bsf hisf endian idf
TESTFILES= getpartest stropentest extstrtest commandtest \
           testio testfs testprompt memiotest mstropentest

help:
	@echo NEMO/src/kernel/io
	@echo bin: BINFILES=$(BINFILES)
	@echo tests: TESTFILES=$(TESTFILES)

clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

#	update or install a library
update:  .update_inc .update_lib 

.update_lib: $(L)
	@touch .update_lib

.update_inc: .install_inc

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
bin: $(BINFILES)
obj: $(OBJFILES)
tests: $(TESTFILES)

hisf: history.c
	$(CC) $(CFLAGS) -DTOOLBOX -o hisf history.c $(NEMO_LIBS)

idf: idfio.c
	$(CC) $(CFLAGS) -DTOOLBOX -o idf idfio.c $(NEMO_LIBS)

endian: endian.c
	$(CC) $(CFLAGS) -DTOOLBOX -o endian endian.c $(NEMO_LIBS)

getpartest: getparam.c
	$(CC) $(CFLAGS) -o getpartest -DTESTBED getparam.c $(NEMO_LIBS)

stropentest: stropen.c
	$(CC) $(CFLAGS) -o stropentest -DTESTBED stropen.c $(NEMO_LIBS)

mstropentest: mstropen.c
	$(CC) $(CFLAGS) -o mstropentest -DTESTBED mstropen.c $(NEMO_LIBS)

extstrtest: extstring.c
	$(CC) $(CFLAGS) -o extstrtest -DTESTBED extstring.c $(NEMO_LIBS)

storytest: story.[ch] storyio.h snode.h
	$(CC) $(CFLAGS) -o storytest -DTOOLBOX story.c $(NEMO_LIBS)

memiotest: memio.c
	$(CC) $(CFLAGS) -o memiotest -DTOOLBOX memio.c $(NEMO_LIBS)

commandtest: command.c
	$(CC) $(CFLAGS) -o commandtest -DTESTBED command.c $(NEMO_LIBS)

# peculiar tests

testio:
	$(CC) $(CFLAGS) -o testio test/testio.c $(NEMO_LIBS)

testi:  getparam.c testi.c
	$(CC) $(CFLAGS) -o testi testi.c $(NEMO_LIBS)

testfs:
	$(CC) $(CFLAGS) -o testfs test/testfs.c $(NEMO_LIBS)

testprompt:
	$(CC) $(CFLAGS) -o testprompt test/testprompt.c $(NEMO_LIBS)

# very special: need root permission to do this

magic:
	cat magic.nemo >> /etc/magic
