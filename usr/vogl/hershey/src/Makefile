#
# Makes the hershey library.
# Makes the binary Hershey font file for VOGL
#
CC = cc
CFLAGS = -I../../src $(MCFLAGS) -DFONTLIB=\"$(FONTLIB)\"

OURLIBS = libhershey.a ../../src/libvogl.a
LIBS = -lsuntool -lsunwindow -lpixrect 
RANLIB = ranlib
#
# On a real SGI ...
#
#OURLIBS = libhershey.a
#LIBS = -lgl
#RANLIB = echo
SHELL = /bin/sh

.IGNORE:

HDR = h2v.h

all:	libhershey.a h2v hdisp fdisp fonts

libhershey.a: htext.o check.o halloc.o fhtext.o
	ar rcv libhershey.a htext.o check.o halloc.o fhtext.o
	$(RANLIB) libhershey.a


h2v:	h2v.o getchar.o
	$(CC) -o $@ $@.o getchar.o

hdisp:	hdisp.o getchar.o ../../src/libvogl.a libhershey.a
	$(CC) -o $@ $@.o getchar.o $(OURLIBS) $(LIBS)

fdisp:	fdisp.o ../../src/libvogl.a libhershey.a
	$(CC) -o $@ $@.o $(OURLIBS) $(LIBS)

$(FONTLIB):
	if test ! -d $(FONTLIB);\
	then mkdir $(FONTLIB);\
	fi;

fonts:	./h2v $(FONTLIB)
	./h2v ../data/hersh.oc
	./h2v ../data/hersh.or ../fonts/japan.hmp japanese
	mv astrology cursive cyrillic futura.l futura.m gothic.eng \
	gothic.ger gothic.ita greek japanese markers math.low \
	math.upp meteorology music script symbolic times.g \
	times.i times.ib times.r times.rb $(FONTLIB)
	touch fonts

h2v: h2v.h

clean:
	rm -f astrology cursive cyrillic futura.l futura.m gothic.eng \
	gothic.ger gothic.ita greek japanese markers math.low \
	math.upp meteorology music script symbolic times.g \
	times.i times.ib times.r times.rb *.o core fonts libhershey.a

clobber:
	rm -f astrology cursive cyrillic futura.l futura.m gothic.eng \
	gothic.ger gothic.ita greek japanese markers math.low \
	math.upp meteorology music script symbolic times.g \
	times.i times.ib times.r times.rb *.o core h2v hdisp fdisp fonts \
	libhershey.a
