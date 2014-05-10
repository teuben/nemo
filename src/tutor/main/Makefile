#	Makefile for a tutorial area

CFLAGS = -g 
FFLAGS = -g -C -u

#	List of all binaries and sources here

BINS = mainf mainc mainc2 nemomain
SRCS = mainf.f mainc.c mainc2.c nemomain.c Makefile

#	Define an extra SUFFIX for .doc file
.SUFFIXES: .doc

.c.doc: $*
	$* help=t > $*.doc
	@echo "### Normally this $*.doc file would be moved to NEMODOC"

#	First target is a dummy, to prevent that is really makes
#	something when 'make' is called without arguments in this
#	directory

help:
	@echo "Tutorial for Nemo programs"
	@echo ""
	@echo "This directory contains a Makefile, with directives"
	@echo "how programs are to be compiled, and the source code"
	@echo "for a fortran, c, and nemo-c program:"
	@echo ""
	@echo "mainf.f		fortran program"
	@echo "mainc.c		c program"
	@echo "nemomain.c	nemo-c program"

#	A few global targets

all: $(BINS)

clean:
	rm -f $(BINS) main core *.o *.doc mainc2.out

#	Install targets - not useful for the person who 
#	studies the tutorial

install: .install_src

.install_src:
	cp $(SRCS) $(NEMO)/src/tutor/main
	@touch .install_src

#	List of individual targets

mainf:	mainf.f
	$(FC) $(FFLAGS) -o mainf mainf.f

mainc:	mainc.c
	$(CC) $(CFLAGS) -o mainc mainc.c

mainc2:	mainc2.c
	$(CC) $(CFLAGS) -I$(NEMOINC) -I$(NEMOLIB) -o mainc2 mainc2.c $(NEMOLIB)/libnemo.a -lm

nemomain: nemomain.c
	$(CC) $(CFLAGS) -o nemomain nemomain.c $(NEMOLIB)/libnemo.a -lm



