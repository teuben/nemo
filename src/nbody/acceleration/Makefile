#	NEMO/src/nbody/core
include $(NEMOLIB)/makedefs

L = $(NEMOLIB)/libnemo.a 
.PRECIOUS: $(L)


MAN3FILES = 
INCFILES = 
SRCFILES = 
OBJFILES = acceleration.o
LOBJFILES = $L(acceleration.o)
BINFILES = 
TESTFILES =

SRCDIR = $(NEMO)/src/nbody/acceleration

help:
	@echo NEMO V2.2 NEMO/src/nbody/acceleration

#install:  .install_lib
install:  
	@echo This directory is deprecated. See src/orbit/potential/acceleration.cc

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	$(RANLIB) $(L)

.install_inc: $(INCFILES)
	cp $? $(NEMO)/inc/snapshot

# Targets used by Nemo to update include libraries and bin directories.

nemo_lib: $(INCFILES)
	@echo no nemo_lib in `pwd`

nemo_bin: $(BINFILES)
	@echo no nemo_bin in `pwd`

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:

clean:
	rm -f *.o *.a *~ core $(BINFILES) $(TESTFILES)

cleanlib:
	ar dv $(L) $(OBJFILES)
tidy:
	rm -f *.o

# special targets

bodytrans:	bodytrans.c
	$(CC) $(CFLAGS) -o bodytrans bodytrans.c -DTOOLBOX \
		$(NEMO_LIBS) $(EL) -lm

# test routines

all: $(BINFILES)  $(TESTFILES)

bin: $(BINFILES)
tests: $(TESTFILES)

