# NEMO 2.x <Makefile>
include $(NEMOLIB)/makedefs

# from: ~josh/hackcode/source1/Makefile: mantain C version of treecode.
#	-o not allowed for .o files ; cleanup up the logic
#	libT added
#	clean *.o afterwards
#
# Notes for speedup:
#  On a Sun:  '-fsingle -DSINGLEPREC -O' speeds up by perhaps 15% AND
#		also takes up 2* less diskspace
#

#CFLAGS = -O
#LDFLAGS = -O
#	Use of -fsingle on sun also required to get real speedup
#CFLAGS = -DSINGLEPREC -O4 -fsingle
#LDFLAGS = -DSINGLEPREC -O4 -fsingle

L= $(NEMOLIB)/libnemo.a
LIBN = $(NEMO_LIBS)
# Files exported to Nemo.

MANFILES = hackcode1.1

BINFILES = hackcode1 hackcode1_qp hackforce_qp hackforce 

SRCFILES = code.c code.h code_io.c defs.h grav.c load.c util.c hackforce.c

SRCDIR = $(NEMOPATH)/src/nbody/evolve/hackcode/source1

# Targets used to export files to Nemo.

install: .install_bin
#install:

.install_man: $(MANFILES)
	cp $? $(NEMOPATH)/man/man1
	@touch .install_man

.install_bin: $(BINFILES)
	cp $? $(NEMOBIN)

.install_src: $(SRCFILES) Makefile
	@if [ ! -d $(SRCDIR) ]; \
	then \
		mkdir $(SRCDIR); \
		chmod 777 $(SRCDIR) ; \
	fi
	cp $? $(SRCDIR)
	@touch .install_src

test:	$(BINFILES)

# Targets used by Nemo to mantain bin files.

nemo_lib:
	@echo no nemo_lib here

nemo_bin: $(BINFILES)
	mv $? $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 
		
		
clean:
	rm -f *.o *.a core $(BINFILES)

tidy:
	rm -f *.o $(BINFILES)
#
hackcode1: code.o code_io.o load_mp.o grav_mp.o util.o
	$(CC) $(CFLAGS) -o hackcode1 \
	   code.o code_io.o load_mp.o grav_mp.o util.o $(LIBN) -lm

hackcode1_qp: code.o code_io.o load_qp.o grav_qp.o util.o
	$(CC) $(CFLAGS) -o hackcode1_qp \
	   code.o code_io.o load_qp.o grav_qp.o util.o $(LIBN) -lm

hackforce: hackforce.c defs.h load_mp.o grav_mp.o
	$(CC) $(CFLAGS) -o hackforce hackforce.c load_mp.o grav_mp.o $(LIBN) -lm

hackforce_qp: hackforce.c defs.h load_qp.o grav_qp.o
	$(CC) $(CFLAGS) -DQUADPOLE -o hackforce_qp \
	   hackforce.c load_qp.o grav_qp.o $(LIBN) -lm

code.o: code.c defs.h code.h

code_io.o: code_io.c defs.h

util.o: util.c defs.h

#	separate out monopole and quadruple versions .o files

load_mp.o: load.c defs.h
	$(CC) $(CFLAGS) -c load.c
	mv load.o load_mp.o

load_qp.o: load.c defs.h
	$(CC) $(CFLAGS) -DQUADPOLE -c load.c
	mv load.o load_qp.o

grav_mp.o: grav.c defs.h
	$(CC) $(CFLAGS) -c grav.c
	mv grav.o grav_mp.o

grav_qp.o: grav.c defs.h
	$(CC) $(CFLAGS) -DQUADPOLE -c grav.c
	mv grav.o grav_qp.o

util.o: util.c defs.h
#
# 2_D versions (not exported to NEMO though)

hackcode1_2d: code_2d.o code_io_2d.o load_2d.o grav_2d.o util_2d.o
	$(CC) $(CFLAGS) -o hackcode1_2d \
	   code_2d.o code_io_2d.o load_2d.o grav_2d.o util_2d.o \
	   $(LIBN) -lm

code_2d.o: code.c defs.h code.h
	$(CC) $(CFLAGS) -DTWODIM -c code.c
	mv code.o code_2d.o

code_io_2d.o: code_io.c defs.h
	$(CC) $(CFLAGS) -DTWODIM -c code_io.c
	mv code_io.o code_io_2d.o

load_2d.o: load.c defs.h
	$(CC) $(CFLAGS) -DTWODIM -DFCELLS=1.0 -c load.c
	mv load.o load_2d.o

grav_2d.o: grav.c defs.h
	$(CC) $(CFLAGS) -DTWODIM -c grav.c
	mv grav.o grav_2d.o

util_2d.o: util.c defs.h
	$(CC) $(CFLAGS) -DTWODIM -c util.c
	mv util.o util_2d.o

hackforce_2d: hackforce.c defs.h load_2d.o grav_2d.o
	$(CC) $(CFLAGS) -DTWODIM -o hackforce_2d hackforce.c load_2d.o grav_2d.o \
	   $(LIBN) -lm

# Misc. test versions.

hackforce_hex: hackforce.c defs.h load.o grav_hex.o
	$(CC) $(CFLAGS) -o hackforce_hex hackforce.c load.o grav_hex.o \
	   $(LIBN) -lm

grav_hex.o: grav_hex.c defs.h
