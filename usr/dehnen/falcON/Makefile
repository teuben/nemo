# -*- makefile -*-
################################################################################
#
# Makefile for the falcON project
#
################################################################################
ifdef NEMO
install		: TestGrav gyrfalcON
else
install		: TestGrav
endif

COMPILER	:= gcc
#COMPILER	:= pgcc
#COMPILER	:= icc

DSOFT 		:= -DfalcON_INDI #allow for individual softening length         

makepublic 	:= Makefile make.defs make.public

include make.defs
include make.public

################################################################################

#pjt hack
test1:	
	g++ -o i386_linux/mkcold src/mains/mkcold.cc -Iinc/ -I/astromake/opt/nemo/cvs/inc -I/astromake/opt/nemo/cvs/lib/i686-pc-linux-gnu -D_FILE_OFFSET_BITS=64 -O3 -ffast-math -funroll-loops -fforce-addr -rdynamic     -DfalcON_NEMO   -Li386_linux/lib/ -lfalcON -lm

nemo_bins:
	make gyrfalcON
	mv $(HOS)/gyrfalcON $(NEMOBIN)
	mv $(LIB)/lib*.so $(NEMOLIB)
