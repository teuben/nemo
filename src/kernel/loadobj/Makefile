#    $NEMO/src/nemo/kernel/loadobj/Makefile
include $(NEMOLIB)/makedefs

# machine name, must have corresponding loadobj<MACH>.c 
# a #ifdef switch table is then present in loadobj.c
# some architectures are automatically loaded, unless -DAUTO is used in CFLAGS

# MACH = 
# MacOS-X seems to need
#    MACH=-DDL -I$(NEMO)/host/$(NEMOHOST)/include
#    EL = -ldl
# extra libraries needed in linking [-lld for SYSV, -ldl for dlopen()]
EL = 

#
# The -D$(MACH) is mainly for loadobj.c, to include the proper loadobj<MACH>.c 
# this needs to fixed up with the makedefs
#CFLAGS = -O -D$(MACH)
#LDFLAGS = -O -D$(MACH)
#RANLIB = ranlib

#	NEMO Library, and make it a precious one
L = $(NEMOLIB)/libnemo.a
.PRECIOUS: $(L)

# Files exported to Nemo.

MAN3FILES = 
INCFILES = loadobj.h
SRCFILES = loadobj.c mapsys.c getfunc.c loadobj*.c
BINFILES = 
OBJFILES = loadobj.o mapsys.o getfunc.o
LOBJFILES = $L(loadobj.o) $L(mapsys.o) $L(getfunc.o)
SRCDIR = $(NEMO)/src/nemo/kernel/loadobj
BINFILES = getfunc
TESTFILES= getfunc loadobjtest mapsystest


# Subdirectories with Nemo Makefiles. 

SUBDIRS = 

# Targets used to export things to Nemo.

help:
	@echo NEMO/src/kernel/loadobj:
	@echo "Targets:"
	@echo "   loadobj, mapsys, getfunc"
	@echo ""
	@echo "To compile a test version of loadobj for a specific machine:"
	@echo "   make loadobjtest MACH=xxx EL=yyy"
	@echo "MACH denotes the loadobj type, as in loadobjMACH.c"
	@echo "EL are extra libraries needed if that MACH type is choosen"
	@echo "MACH=$(MACH)   EL=$(EL)"

clean:
	rm -f *.o *.so *.a *~ core $(BINFILES) $(TESTFILES) loadobj mapsys ld-tmp.* ms-tmp.*

cleanlib:
	ar dv $(L) $(OBJFILES)

all:	$(BINFILES) $(TESTFILES)

# Targets used to assemble object library.

tests: $(TESTFILES)

#
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

bins:	$(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

####    Test with DLOPEN()  - sep-90 - pjt

loadobjtest: loadobj.c testloadobj.c 
	$(CC) $(CFLAGS) -o loadobjtest -DTESTBED \
		loadobj.c $(NEMO_LIBS) $(EL) 

getfunc: getfunc.c
	$(CC) $(CFLAGS) -o getfunc -DTOOLBOX \
		getfunc.c $(NEMO_LIBS) $(EL) 

mapsystest: mapsys.c
	$(CC) $(CFLAGS) -o mapsystest -DTESTBED \
		mapsys.c $(NEMO_LIBS) $(EL) 

#	dependencies on loadobj.o are not quite as they should be
#	since loadobj.c depends on all the loadobj*.c files
#
# loadobj.o: loadobj.c loadobj.h
#	$(CC) $(CFLAGS) -D$(MACH) -c loadobj.c

#	special make entries

DLD=dld-3.2.3

dld:
	(cd $(NEMO)/usr/lib/$(DLD); \
	 $(MAKE) libdld.a CFLAGS="$(CFLAGS) CC="$(CC)" ;\
	 mv libdld.a $(NEMOLIB) ; \
         $(MAKE) clean)
