#
# examples makefile 
#
CC = cc

CEXAMPS = trivial simple shapes poly views circtxt moretxt curves \
	  moretxt2 patches balls objvws world loc teapot cube lcube tetra \
	  sinewave tmesh

COBJS = moretxt2.o trivial.o simple.o shapes.o poly.o views.o circtxt.o moretxt.o \
	curves.o patches.o balls.o objvws.o world.o loc.o \
	teapot.o cube.o lcube.o tetra.o sinewave.o tmesh.o

#
# Where to find librarys
#
OURLIBS = ../src/libvogl.a ../hershey/src/libhershey.a
LIBS = -lX11 -lsuntool -lsunwindow -lpixrect -lm
MCFLAGS = -g -fsingle /usr/lib/libm.il
#
# On a real SGI...
#
#OURLIBS = ../hershey/src/libhershey.a
#LIBS = -lgl
#MCFLAGS = -g -float

CFLAGS = -I../src $(MCFLAGS)

all:	$(CEXAMPS)

$(CEXAMPS): $(COBJS) $(OURLIBS)
	$(CC) $(CFLAGS) -o $@ $@.o $(OURLIBS) $(LIBS)

clean:
	rm -f *.o core

clobber:
	rm -f $(CEXAMPS) *.o core
