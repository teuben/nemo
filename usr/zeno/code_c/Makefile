##############################################################################
# MAKEFILE: compile tree-code program.                                       #
#                                                                            #
# Copyright (c) 1993, Joshua E. Barnes, Honolulu, HI.                        #
# It's free because it's yours.                                              #
##############################################################################

# Compiler flags for development on SGI; not used on other systems.
# single version: CFLAGS = -fullwarn -woff 269,302 -prototypes -O -float
# double version: CFLAGS = -DDOUBLEPREC -fullwarn -woff 269,302 -prototypes -O

# Standard compiler flags.  Double precision is the default because not
# all available systems seem to have single-precision square root.
CFLAGS = -DDOUBLEPREC -O

all:   code

clean:
	rm -f *.o code

# Build a treecode.

code: code.o io.o load.o grav.o util.o getparam.o
	$(CC) $(LDFLAGS) -o code code.o io.o load.o grav.o\
		util.o getparam.o -lm

code.o: code.c defs.h code.h stdinc.h
	$(CC) $(CFLAGS) -c code.c

io.o: io.c defs.h stdinc.h
	$(CC) $(CFLAGS) -c io.c

load.o: load.c defs.h stdinc.h
	$(CC) $(CFLAGS) -c load.c

grav.o: grav.c defs.h stdinc.h
	$(CC) $(CFLAGS) -c grav.c

#  Utility routines etc.

util.o: util.c stdinc.h
	$(CC) $(CFLAGS) -c util.c

getparam.o: getparam.c stdinc.h
	$(CC) $(CFLAGS) -c getparam.c

