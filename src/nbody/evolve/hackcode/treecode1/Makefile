##############################################################################
# MAKEFILE: compile and link treecode executable files                       #
# Copyright (c) 2001 by Joshua E. Barnes, Honolulu, Hawai`i.                 #
##############################################################################

########################################################################
# DIRECTIONS: to tune this Makefile to your system, edit the definitions
# of CCFLAGS, LDFLAGS, OPTFLAG, and PREC below.  Suggested values
# are provided for LINUX, Sun, and SGI systems.

########################################################################
# Compiler options.

# LINUX:
CCFLAGS = -DLINUX
LDFLAGS =
OPTFLAG = -O3

# Sun:
# CCFLAGS =
# LDFLAGS =
# OPTFLAG = -xO3

# SGI:
# CCFLAGS = -n32
# LDFLAGS = -n32
# OPTFLAG = -O3

########################################################################
# Precision.  Possible values are SINGLEPREC, MIXEDPREC, and DOUBLEPREC.

# LINUX, SGI:
PREC = SINGLEPREC

# Sun:
# PREC = MIXEDPREC

########################################################################
# Variations:

# Uncomment the next line to use freq instead of dtime:
# USEFREQ = -DUSEFREQ

# Uncomment the next line to use binary I/O:
# BINARYIO = -DBINARYIO

########################################################################
# Vanila treecode:

treecode: treecode.o treeio.o treeload.o treegrav.o libZeno.a
	$(CC) $(LDFLAGS) -o treecode \
	  treecode.o treeio.o treeload.o treegrav.o libZeno.a -lm

treecode.o: treecode.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) $(USEFREQ) -c treecode.c

treeio.o: treeio.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) $(USEFREQ) $(BINARYIO) -c treeio.c

treeload.o: treeload.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) -c treeload.c

treegrav.o: treegrav.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) $(OPTFLAG) -c treegrav.c

########################################################################
# Quick-scan treecode:

treecode_q: treecode_q.o treeio_q.o treeload_q.o treegrav_q.o libZeno.a
	$(CC) $(LDFLAGS) -o treecode_q \
	  treecode_q.o treeio_q.o treeload_q.o treegrav_q.o libZeno.a -lm

treecode_q.o: treecode.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) -DQUICKSCAN $(USEFREQ) -c \
	  -o treecode_q.o treecode.c

treeio_q.o: treeio.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) -DQUICKSCAN $(USEFREQ) $(BINARYIO) -c \
	  -o treeio_q.o treeio.c

treeload_q.o: treeload.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) -DQUICKSCAN -c \
	  -o treeload_q.o treeload.c

treegrav_q.o: treegrav.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PREC) -DQUICKSCAN $(OPTFLAG) -c \
	  -o treegrav_q.o treegrav.c

########################################################################
# Zeno library:

libZeno.a: clib.o getparam.o mathfns.o
	ar ruv libZeno.a clib.o getparam.o mathfns.o

clib.o: clib.c stdinc.h
	$(CC) $(CCFLAGS) -D$(PREC) -c clib.c

getparam.o: getparam.c stdinc.h getparam.h
	$(CC) $(CCFLAGS) -D$(PREC) -c getparam.c

mathfns.o: mathfns.c stdinc.h mathfns.h
	$(CC) $(CCFLAGS) -D$(PREC) -c mathfns.c
