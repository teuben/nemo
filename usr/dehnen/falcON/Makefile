# -*- makefile -*-
################################################################################
#
# Makefile for the falcON project
#
################################################################################
default		: all
ifdef NEMO
tonemo		: all
			-cp $(LIB)/libfalcON.so $(NEMOLIB)
			-cp $(BIN)/* $(NEMOBIN) >& /dev/null
install		: tonemo man
else
install		: all man
endif

COMPILER	:= gcc
#COMPILER	:= pgcc
#COMPILER	:= icc

DSOFT 		:= -DfalcON_INDI #allow for individual softening length         

makepublic 	:= Makefile make.defs make.public

include make.defs
include make.public

################################################################################
# pjt (probably not needed anymore now)
ifdef NEMO
include Makefile.nemo
endif
