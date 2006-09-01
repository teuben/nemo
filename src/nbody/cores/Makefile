#	NEMO/src/nbody/core
include $(NEMOLIB)/makedefs

L = $(NEMOLIB)/libnemo.a 
.PRECIOUS: $(L)


MAN3FILES = 
INCFILES = snapshot.h \
	   barebody.h body.h mybody.h \
	   stdbody.h \
	   units.h
SRCFILES = snapshot.h barebody.h body.h get_snap.c put_snap.c snaptest.c
OBJFILES = pickpnt.o units.o zerocms.o bodytrans.o
LOBJFILES = $L(pickpnt.o) $L(units.o) $L(zerocms.o) $L(bodytrans.o)
BINFILES = bodytrans
TESTFILES = testunits

SRCDIR = $(NEMO)/src/nbody/io

help:
	@echo NEMO V2.2 NEMO/src/nbody/cores

install:  .install_lib

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

testunits: units.[ch]
	$(CC) $(CFLAGS) -o testunits -DTESTBED units.c $(NEMO_LIBS)
