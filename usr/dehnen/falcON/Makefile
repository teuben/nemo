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
