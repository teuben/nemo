########################################################################
# MAKEFILE: compile and link treesph code executable files             #
# Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                #
#    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000       #
########################################################################


############################
### Compiler and Options ###
############################

# Linux, Sun:
#CC = fcc
CC = gcc
CCFLAGS = -O3
LDFLAGS = 

# Digital UNIX:
# CC = cc
# CCFLAGS = -O3
# LDFLAGS =

#########################################################
### Precision. [SINGLEPREC, MIXEDPREC, or DOUBLEPREC] ###
#########################################################

# Linux, Sun, Digital UNIX:
PRECISION = MIXEDPREC

#############
### Links ###
#############

X11_LIB = -L/usr/X11R6/lib -lX11

F2C_LIB = -lg2c

CPGPLOT_LIB = -L/usr/local/pgplot/ -lcpgplot -lpgplot -lpng -lgcc -lz\
         $(X11_LIB) $(F2C_LIB)

#####################
### Tree-Sph code ###
#####################


OBJECTS = treecode.o treeio.o treeload.o treegrav.o treengbr.o \
	treesph.o treeextf.o treeinit.o libZeno.a

etude: $(OBJECTS) libZeno.a
	$(CC) $(LDFLAGS) -o etude $(OBJECTS) -lm

treecode.o: treecode.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treecode.c

treeio.o: treeio.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -DBINARYIO -c treeio.c

treeload.o: treeload.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treeload.c

treegrav.o: treegrav.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treegrav.c

treengbr.o: treengbr.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treengbr.c

treesph.o: treesph.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treesph.c

treeextf.o: treeextf.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treeextf.c

treeinit.o: treeinit.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c treeinit.c


############################################
### Tree-Sph code with quick look option ###
############################################

OBJECT_QL = treecode_ql.o treeio.o treeload.o treegrav.o treengbr.o \
		treesph.o treeextf.o treeinit.o plot.o

etude_ql: $(OBJECT_QL) libZeno.a
	$(CC) $(LDFLAGS) -o etude_ql $(OBJECT_QL) libZeno.a \
	-lm $(CPGPLOT_LIB)

treecode_ql.o: treecode.c treecode.h treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -DQLOOK -c \
	-o treecode_ql.o treecode.c 

plot.o: plot.c treedefs.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c plot.c

####################
### Zeno library ###
####################

libZeno.a: clib.o getparam.o mathfns.o
	ar ruv libZeno.a clib.o getparam.o mathfns.o
	ranlib libZeno.a

clib.o: clib.c stdinc.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c clib.c

getparam.o: getparam.c stdinc.h getparam.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c getparam.c

mathfns.o: mathfns.c stdinc.h mathfns.h
	$(CC) $(CCFLAGS) -D$(PRECISION) -c mathfns.c

#####################
### Data analysis ###
#####################

movie: movie.c stdinc.h mathfns.h vectmath.h
	$(CC) -o movie movie.c  -lm $(CPGPLOT_LIB)

analysis: analysis.c stdinc.h mathfns.h vectmath.h
	$(CC) -o analysis analysis.c -lm $(CPGPLOT_LIB)

###########################
### Remove object files ###
###########################

clean:;	rm -f *.o libZeno.a
