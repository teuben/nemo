#
# examples makefile 
#
FEXAMPS = ftrivial fsimple fshapes fpoly fviews fcirctxt fmoretxt fmoretx2 \
	  fcurves fpatches fballs fobjvws fworld floc ftetra fcube \
	  flcube fsinwave ftmesh

FOBJS = ftrivial.o fsimple.o fshapes.o fpoly.o fviews.o fcirctxt.o fmoretxt.o \
	fmoretx2.o fcurves.o fpatches.o fballs.o fobjvws.o fworld.o floc.o \
	ftetra.o fcube.o flcube.o fsinwave.o ftmesh.o

#
# Where to find library
LIB = ../hershey/src/libhershey.a ../src/libvogl.a
LIBS = -lsuntool -lsunwindow -lpixrect
#

.SUFFIXES: .F .o

F77 = f77

MFFLAGS = -g /usr/lib/libm.il
FFLAGS = $(MFFLAGS) -I../src

.F.o:
	$(F77) -c $(FFLAGS) $*.F

all:	$(FEXAMPS)

$(FEXAMPS): $(FOBJS) $(LIB)
	$(F77) $(FFLAGS) -o $@ $@.o $(LIB) $(LIBS) 

clean:
	rm -f *.o core

clobber:
	rm -f $(FEXAMPS) *.o core
