#
# Makefile for vogl Set up for SUN Solaris2.x
#
# Usage: make
#
# As well as selecting the devices, you may have to change LLIBS and MLIBS 
# as well (see below).
#
# Some machines use /bin/csh in Make (or use the envionment variable
# SHELL. We want to us /bin/sh.
#
SHELL = /bin/sh
#
#  Which cc to use
#
#CC = gcc
CC = cc

# 
#  The devices you wish to have compiled into the library and the
#  location of the object files for each device relative to the 
#  src directory. For each device defined there should be a 
#  corresponding ../drivers/file.o in the DOBJS line.
#
#  Possible DEVICES and their object files are:
#		-DPOSTSCRIPT	../drivers/ps.o
#		-DHPGL		../drivers/hpdxy.o
#		-DDXY		../drivers/hpdxy.o
#		-DTEK		../drivers/tek.o
#		-DSUN		../drivers/sun.o	(Sunview that is)
#		-DX11		../drivers/X11.o
#		-DNeXT		../drivers/NeXT.o	(NeXTStep)
#
DEVICES = -DPOSTSCRIPT -DHPGL -DTEK -DX11
DOBJS = ../drivers/ps.o ../drivers/hpdxy.o ../drivers/tek.o ../drivers/X11.o
#
# Where the fonts a going to live (For making the Hershey library)
# (Remember /tmp usually gets cleared from time to time)
#
FONTLIB = /usr/local/lib/hershey/
FONTLIB = /tmp/hfonts

#
#  For BSD type machines we need to use ranlib
#
#RANLIB = ranlib
# If you have SYS5 then set RANLIB = "ar ts" or maybe "echo" in later versions
RANLIB = echo

#
# Set any Special floating point options here...
# 
# The default 
FLOATING_POINT =

#
# On a SUN (use inline math functions, don't promote to double in expressions)
FLOATING_POINT = -xlibmil -fsingle

# On a SUN3 with an mc68881 floating point chip
#FLOATING_POINT = /usr/lib/libm.il -fsingle -f68881

#
# On a NeXT 	(for inline expansion of math primitives)
#FLOATING_POINT = -DINLINE_MATH

# On an SGI
#FLOATING_POINT = -float

# Or on a Decstation (running Ultrix of course)
#FLOATING_POINT = -f

#
# Global CFLAGS can be set here.
#
# The default
#CFLAGS = -O -I/usr/local/R5/include

#
#  SUN4/Solaris2
CFLAGS = -DSYS5 -xO4 -I/usr/openwin/include

# For SUNOS 3.5 set -DSUN_3_5
#
#CFLAGS = -O4 -DSUN_3_5

#
# For SGI's
#
#CFLAGS = -DSYS5 -O

#
# Or on an apollo... (-Wp make it use the 'other cpp')
#CFLAGS = -O -Wp -M3000

#
# Or an IBM RS6000
#
#CFLAGS = $(OS) -O -Q

#
# NeXTStep....
#
#CFLAGS = -O -ObjC

#
# Define F77 if you want the f77 examples.
F77 = f77
#
# You also define your f77 flags here too. These are the ones we used on sun
#
FFLAGS = -O -w -libmil
#
# Or on an apollo (We didn't have ftn)
#F77 = 
#FFLAGS = 

#
# Or on an IBM RS6000
#
#F77 = xlf
#FFLAGS = -O -qextname

# The usual default...
#FFLAGS = -O


#
# The name of the library to install and where to put it.
#
LIB = libvogl.a
DEST = /usr/local/lib


#
# Any other local libraries that are needed go here
#

# SUN/Solaris
LLIBS = -lX11 -L/usr/openwin/lib
# X11
#LLIBS = -lX11
# SGI
#LLIBS = -lgl_s
# NeXT
#LLIBS = -lNeXT_s

#
# On the NeXT there is some kind of problem with using -lm when
# linking. Leave LIBM blank in this case.
#

LIBM = -lm

LIBS = $(LLIBS) $(LIBM)

MCFLAGS = $(CFLAGS) $(FLOATING_POINT)
MFFLAGS = $(FFLAGS)

all:
	cd src; make -f Makefile \
			CC="$(CC)" \
			DEVICES="$(DEVICES)" \
			MCFLAGS="$(MCFLAGS)" \
			DOBJS="$(DOBJS)"\
			RANLIB="$(RANLIB)"

	cd hershey/src; make -f Makefile \
			CC="$(CC)" \
			FONTLIB="$(FONTLIB)" \
			MCFLAGS="$(MCFLAGS)" \
			LIBS="$(LIBS)" \
			RANLIB="$(RANLIB)"

	cd examples; make -f Makefile \
			CC="$(CC)" \
			MCFLAGS="$(MCFLAGS)" \
			LIBS="$(LIBS)"

	cd examples/xview; make -f Makefile \
			CC="$(CC)" \
			MCFLAGS="$(MCFLAGS)" \
			LIBS="$(LIBS)"

	if test -n "$(F77)" ; \
	then cd examples; make -f Makefile.f77 \
			LIBS="$(LIBS)" \
			MFFLAGS="$(MFFLAGS)" \
			F77="$(F77)" ; \
	fi ; exit 0

install:
	cp src/$(LIB) $(DEST)
	chmod 644 $(DEST)/$(LIB)
	$(RANLIB) $(DEST)/$(LIB)

clean:
	cd src; make DOBJS="$(DOBJS)" clean
	cd hershey/src; make FONTLIB="$(FONTLIB)" clean
	cd drivers; make clean
	cd examples; make clean
	cd examples/xt; make clean
	cd examples/xview; make clean
	cd examples/sunview; make clean
	cd examples; make -f Makefile.f77 clean

clobber:
	cd src; make DOBJS="$(DOBJS)" clobber
	cd hershey/src; make FONTLIB="$(FONTLIB)" clobber
	cd drivers; make clobber
	cd examples; make clobber
	cd examples/xt; make clobber
	cd examples/xview; make clobber
	cd examples/sunview; make clobber
	cd examples; make -f Makefile.f77 clobber
