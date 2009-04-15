# ~josh/multicode/Makefile: maintain multipole codes.
include $(NEMOLIB)/makedefs
#	22-nov-89	HOSTTYPE
#	12-nov-91	new Nemo V2. : location now in new $NEMO/src
#			also deleted code_io.c and defs.h
MBODY=100000

QFLAGS =  -DMBODY=$(MBODY)

# new lib: (YAPPLIB is imported from NEMORC file)
L = $(NEMO_LIBS)
# old lib:
#LIBN = $(NEMOLIB)/libNbody.a
#LIBJ = $(NEMOLIB)/libT.a $(NEMOLIB)/libJ.a 
#LIBG = -lcore68881 -lsunwindow -lpixrect

# Objects exported to NEMO

MAN1FILES = quadcode.1

INCFILES = quadfield.h quadforce.c orbstep.c

BINFILES = quadcode quadforce quadinter

SRCFILES = quaddefs.h quadfield.h quadcode.c quadcode_io.c orbstep.c \
	   quadforce.c quadforce_main.c quadinter.c quadinter_main.c

# SRCDIR = $(NEMO)/src/josh/multicode
SRCDIR = $(NEMO)/src/nbody/evolve/multicode

help:
	@echo Multicode: $(BINFILES)
	@echo MBODY=$(MBODY)

# Targets used to install multipole codes in Nemo.

#install: .install_inc .install_bin .install_src
install:

.install_man: $(MAN1FILES)
	cp $(MAN1FILES) $(NEMO)/man/man1
	@touch .install_man

.install_inc: $(INCFILES)
	cp $(INCFILES) $(NEMO)/inc/multicode
	@touch .install_inc

.install_bin: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	@touch .install_bin

.install_src: $(SRCFILES) $(INCFILES) Makefile
	@if [ ! -d $(SRCDIR) ]; \
	then \
		mkdir $(SRCDIR); \
		chmod 777 $(SRCDIR); \
	fi
	cp $? $(SRCDIR)
	@touch .install_src

# Targets used by Nemo to update include libraries and bin directories.

nemo_lib: 
	@echo No nemo_lib here

nemo_bin: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

clean:
	rm -f *.o *.a *~ core $(BINFILES)
#
tidy:
	rm -f *.o *.BAK *.CKP

all:	$(BINFILES)

quadcode: quadcode.o orbstep.o quadforce.o quadcode_io.o
	$(CC) $(CFLAGS)  $(QFLAGS) -o quadcode quadcode.o orbstep.o quadforce.o \
	    quadcode_io.o $(L) -lm

quadcode.o: quadcode.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadcode.c

orbstep.o: orbstep.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c orbstep.c

quadforce.o: quadforce.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadforce.c

quadcode_io.o: quadcode_io.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadcode_io.c

quadforce: quadforce_main.o quadforce.o
	$(CC) $(CFLAGS)  $(QFLAGS) -o quadforce quadforce_main.o quadforce.o \
	    $(L) -lm

quadforce_main.o: quadforce_main.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadforce_main.c

quadinter: quadinter_main.o quadinter.o
	$(CC) $(CFLAGS) $(QFLAGS)  -o quadinter quadinter_main.o quadinter.o \
	    $(L) -lm

quadinter_main.o: quadinter_main.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadinter_main.c

quadinter.o: quadinter.c quaddefs.h quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -c quadinter.c

# Misc. test routines. (JEB: development)

orbtest1: orbtest1.c orbstep.o quaddefs.h
	$(CC) $(CFLAGS)  $(QFLAGS) -o orbtest1 orbtest1.c orbstep.o $(L) -lm

orbtest2: orbtest2.c orbstep.o quaddefs.h
	$(CC) $(CFLAGS) $(QFLAGS)  -o orbtest2 orbtest2.c orbstep.o $(L) -lm

q22plot: q22plot.c quadfield.h
	$(CC) $(CFLAGS)  $(QFLAGS) -o q22plot q22plot.c $(L) $(YAPPLIB) -lm
