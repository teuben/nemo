#	NEMO/src/nbody/io
include $(NEMOLIB)/makedefs

SUBDIRS = starlab

L = $(NEMOLIB)/libnemo.a 

# Objects exported to nemo.

MAN3FILES = 
INCFILES = get_snap.c put_snap.c myget_snap.c myput_snap.c \
	   get_snapshot.c put_snapshot.c \
           barebody.h body.h mybody.h snapshot.h sphbody.h
SRCFILES = $(INCFILES)
OBJFILES = snapserial.o
SRCDIR = $(NEMO)/src/nbody/io
SUBDIRS= gadget
BINFILES = atos atos_sp atosph atosph_sp stoa stoa_sp tabtos \
	   eaglesnap \
	   rvsnap snaprv snaptipsy tipsysnap \
	   gadgetsnap snapgadget binsnap

TESTFILES = testss testio

test:
	@echo NEMO NEMO/src/nbody/io

archaic:  .archaic_inc

.archaic_inc:
	cp archaic/*.[ch] $(NEMO)/inc/archaic

install: .install_lib

.install_inc: $(INCFILES)
	cp $? $(NEMOINC)/snapshot

.install_lib: $(OBJFILES)
	ar ruv $(L) $?
	ranlib $(L)

update:


# Targets used by Nemo to update include libraries and bin directories.

nemo_lib: $(INCFILES)
	@echo no nemo_lib in `pwd`

nemo_bin: nemo_subbins $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES) nemo_subbins
	mv $(BINFILES) $(NEMOBIN)

nemo_subbins:
	-@for i in $(SUBDIRS); do \
	(cd $$i; echo SUBSUBDIR $$i; $(MAKE) nemo_bin); done

nemo_src:

clean:
	-@for i in $(SUBDIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) clean); done
	rm -f *.o *.a core $(BINFILES) $(TESTFILES)

all: $(BINFILES)

bin: $(BINFILES)
test: $(TESTFILES)

# Test routine for get_snap, put_snap.

testss: testss.c snapshot.h body.h myget_snap.c myput_snap.c
	$(CC) $(CFLAGS) -o testss testss.c $(NEMO_LIBS)

# warning: old target, not portable
testio:	get_snapshot.c testio.c bodywork.o
	$(CC) $(CFLAGS) -o testio testio.c bodywork.o $(L) \
	  -L/usr/lang/SC0.0 -lF77 -lm

bodywork.o: bodywork.f
	$(FC) $(FFLAGS) -c bodywork.f

#  IO import/export programs

atos: atos.c
	$(CC) $(CFLAGS) -o atos atos.c $(NEMO_LIBS)

atos_sp: atos.c
	$(CC) $(CFLAGS) -o atos_sp -DSINGLEPREC atos.c $(NEMO_LIBS)

atosph: atos.c
	$(CC) $(CFLAGS) -o atosph -DSPH atos.c $(NEMO_LIBS)

atosph_sp: atos.c
	$(CC) $(CFLAGS) -o atosph_sp -DSPH -DSINGLEPREC atos.c $(NEMO_LIBS)

stoa: stoa.c
	$(CC) $(CFLAGS) -o stoa stoa.c  $(NEMO_LIBS)

stoa_sp: stoa.c
	$(CC) $(CFLAGS) -o stoa_sp -DSINGLEPREC stoa.c  $(NEMO_LIBS)

tabtos: tabtos.c
	$(CC) $(CFLAGS) -o tabtos tabtos.c  $(NEMO_LIBS)

xvpsnap: xvpsnap.c
	$(CC) $(CFLAGS) -o xvpsnap xvpsnap.c  $(NEMO_LIBS)

snap3dv: snap3dv.c
	$(CC) $(CFLAGS) -o snap3dv snap3dv.c  $(NEMO_LIBS)

snaptipsy: snaptipsy.c
	$(CC) $(CFLAGS) -o snaptipsy snaptipsy.c  $(NEMO_LIBS)

tipsysnap: tipsysnap.c
	$(CC) $(CFLAGS) -o tipsysnap tipsysnap.c  $(NEMO_LIBS)

# This needs:  (cd $NEMO/usr/eagle; make install)
# @todo  this needs to be done  via configure
#        the current HDF5_LIB and HDF5_INC are not good enough
#H_FLAGS = `pkg-config --cflags hdf5`
#H_LIBS  = `pkg-config --libs hdf5`

eaglesnap: eaglesnap.c
	$(CC) $(CFLAGS) $(HDF5_INC) -o eaglesnap eaglesnap.c -lread_eagle $(HDF5_LIB) $(NEMO_LIBS)
