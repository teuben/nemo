# Makefile for Aarseth NBODY tools in NEMO environment
include $(NEMOLIB)/makedefs

SRCDIR = $(NEMO)/src/nbody/evolve/aarseth/tools
SRC = Makefile nbody_io.f nbody_io_c.c nbody_io.h \
      stou4.c u3tos.c u4tos.c
MAN = stou4.1 u3tos.1 u4tos.1

#
L = $(NEMOLIB)/libnemo.a
OBJFILES= nbody_io.o nbody_io_c.o
BINFILES= u3tos u4tos stou4 runbody1 runbody2 runbody4
TESTFILES= runbody1 runbody2 
#	Define an extra SUFFIX for our .doc file
.SUFFIXES: .doc

.c.doc: $*
	$* help=t > $*.doc
	@echo "### Normally this $*.doc file would be moved to NEMODOC"
	@echo "### You can also use mkpdoc to move it over"

help:
	@echo NEMO V2.2
	@echo "Aarseth NBODY data conversion tools for NEMO"
	@echo " No more help to this date"
	@echo "BINFILES=$(BINFILES)"

#  this rule is not how NEMO normally operates.....

install: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

all:	$(BINFILES) $(TESTFILES)

clean:
	rm -f core *.o *.a *.doc $(BINFILES) $(TESTFILES)

test:
	$(MAKE) -f Testfile all

cleanlib:
	ar dv $(L) $(OBJFILES)
	ranlib $(L)

$(L):   $(LOBJFILES)
	echo "*** Now updating all members ***"
	ar ruv $(L) $?
	$(RANLIB) $(L)
	rm -f $?
				
lib:   $(L)

bin: $(BINFILES)


bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

BL = nbody_io.o

# update from source into NEMO

update: .update_src .update_man

.update_src: $(SRC)
	cp $? $(SRCDIR)
	@touch .update_src

.update_man: $(MAN)
	cp $? $(NEMO)/man/man1
	@touch .update_man

# compiling a C main with FC may not work on all architectures
# in those cases link with CC and find out the fortran libraries
# needed to link

# may also need to switch some FC <-> CC, i.e. it's important who does the linking
# some compilers like to see -nofor_main 

u3tos_f: u3tos.o nbody_io.o
	$(FC) $(F77_CMAIN_FFLAGS) $(FFLAGS) -o u3tos u3tos.o nbody_io.o  $(NEMO_LIBS)

u3tos: u3tos.o nbody_io_c.o
	$(CC) $(CFLAGS) -o u3tos u3tos.o nbody_io_c.o $(NEMO_LIBS)

u4tos: u4tos.o nbody_io.o
	$(FC) $(F77_CMAIN_FFLAGS) $(FFLAGS) -o u4tos u4tos.o nbody_io.o  $(NEMO_LIBS)

stou4: stou4.o nbody_io.o
	$(FC) $(F77_CMAIN_FFLAGS) $(FFLAGS) -o stou4 stou4.o nbody_io.o  $(NEMO_LIBS)

runbody1: runbody1.o nbody_io.o
	$(FC) $(F77_CMAIN_FFLAGS) $(FFLAGS) -o runbody1 runbody1.o nbody_io.o  $(NEMO_LIBS)

runbody2: runbody2.o
	$(CC) $(CFLAGS) -o runbody2 runbody2.o  $(NEMO_LIBS)

