#	Makefile for NEMO/src/nbody/xyz
include $(NEMOLIB)/makedefs
#	(see also Makefile.pjt for some odd machine versions)

#  note that xyzview NEEDS to be in singleprecision (or whatever GL default is)
#  you probably don't want to touch this flag
#  For older VOGL's (pre 1.3.0 ??) you may need an addition -DOLD_VOGL
#  Could perhaps also use -DXMAXSCREEN=512 -DYMAXSCREEN=512
XYZFLAG = -DXYZVEL 

#  for SGI you probably need to use:
#		VOGLINC= VOGLLIB=-lgl
VOGLINC = -I$(NEMOINC)/vogl
VOGLLIB = $(NEMOLIB)/libvogl.a $(XLIBS)

MESA = $(HOME)/code/Mesa-6.0.1
MESAINC = -DMESA -I$(MESA)/include
MESALIB = -L$(MESA)/lib -lGL -lGLU -lglut $(XLIBS)


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
BINFILES = snapxyz xyzview
TESTFILES= 

help:
	@echo NEMO V2.2 NEMO/src/nbody/xyz

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
	echo "*** Now updating all members ***"
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

snapxyz:  snapxyz.c
	$(CC) $(CFLAGS)  $(XYZFLAG) \
		-o snapxyz snapxyz.c \
		$(NEMO_LIBS) $(EL) -lm

xyzview:  xyzview.c
	$(CC) -DSINGLEPREC $(CFLAGS)  $(XYZFLAG)  $(VOGLINC) \
		-o xyzview xyzview.c \
		$(NEMO_LIBS) $(VOGLLIB) $(EL) -lm


mesa:	
	$(CC) -DSINGLEPREC $(CFLAGS)  $(XYZFLAG)  $(MESAINC) \
		-o xyzview xyzview.c \
		$(NEMO_LIBS) $(MESALIB) $(EL) -lm


