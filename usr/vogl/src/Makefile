#
# Make the vogl library.
#
# Below are our defaults for compliers etc. These are set in the
# main Makefile.
#

CC = cc

# Hack because our version of gcc produces code which the gnu as 
# program can't eat. Until someone installs the patcehs then we'll
# use the normal sun as.
#
#.c.o:
#	$(CC) -S $(CFLAGS) $*.c
#	as -o $*.o $*.s

MCFLAGS = -g -fsingle

LIB = libvogl.a

RANLIB = ranlib		# = "ar ts" for SYS5

#DEVICES = -DSUN -DTEK -DHPGL -DPOSTSCRIPT -DPOSTSCRIPTP -DDXY
DEVICES = -DSUN

DOBJS = ../drivers/sun.o

HDRS =	vogl.h

CSRC =	arcs.c attr.c buffer.c clip.c curves.c draw.c drivers.c \
	viewing.c viewp.c mapping.c matrix.c move.c objects.c patches.c \
	points.c polygons.c rect.c tensor.c text.c newtokens.c trans.c \
	scale.c getgp.c valloc.c yobbarays.c pref.c vcalls.c \
	queue.c lines.c

COBJS =	arcs.o attr.o buffer.o clip.o curves.o draw.o drivers.o \
	viewing.o viewp.o mapping.o matrix.o move.o objects.o patches.o \
	points.o polygons.o rect.o tensor.o text.o newtokens.o trans.o \
	scale.o getgp.o valloc.o yobbarays.o pref.o vcalls.o \
	queue.o lines.o


FOBJS =	sunfort/farcs.o sunfort/fattr.o sunfort/fcurves.o \
	sunfort/fdraw.o sunfort/fdrivers.o sunfort/fgetgp.o  \
	sunfort/fviewing.o sunfort/fviewp.o sunfort/fmatrix.o sunfort/fmove.o \
	sunfort/fobjects.o sunfort/fpatches.o sunfort/fpoints.o sunfort/fpoly.o\
	sunfort/frect.o sunfort/ftrans.o sunfort/fyobbarays.o sunfort/ftext.o \
	sunfort/fbuffer.o sunfort/fpref.o sunfort/fqueue.o sunfort/fvcalls.o \
	sunfort/fbgnend.o sunfort/flines.o

CFLAGS= $(DEVICES) $(MCFLAGS)
LINTFLAGS=


$(LIB): $(DOBJS) $(COBJS) $(FOBJS) vogl.h
	rm -f $@
	ar rcv $@ $(DOBJS) $(COBJS) $(FOBJS)
	$(RANLIB) $@

lint:
	lint -hac $(DEVICES) $(CSRC) sunfort/*.c -I. -lm > vogl.lint

clean:
	rm -f $(COBJS) core
	cd sunfort; make clean
	cd ../drivers; make clean

clobber:
	rm -f $(LIB) $(COBJS) $(FOBJS) core *.lint
	cd sunfort; make clean

$(COBJS): vogl.h

$(FOBJS): vogl.h
	cd sunfort; make -f Makefile MCFLAGS="$(MCFLAGS)"

$(DOBJS): vogl.h
	cd ../drivers; make -f Makefile \
			"CC=$(CC)" \
			MCFLAGS="$(MCFLAGS)" \
			DOBJS="$(DOBJS)"

