# -*- makefile -*-
################################################################################
#
# Makefile for the falcON project
#
################################################################################
default		: utils all manips acc doxu
ifdef NEMO
tonemo		: all manips
			cp $(LIB)/libfalcON.so $(NEMOLIB)
		      - cp $(BIN)/* $(NEMOBIN) >& /dev/null
install		: tonemo man
else
install		: default man
endif

COMPILER	:= gcc
#COMPILER	:= pgcc
#COMPILER	:= icc

DSOFT 		:= -DfalcON_INDI# allow for individual softening length         

makepublic 	:= Makefile make/defs make/public make/manip

include make/defs
include make/public
include make/manip
include make/acc

falcON_objs	 =	$(public_objs)

$(LIB)libfalcON.a:	$(falcON_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)libfalcON.so:	$(falcON_objs)
			$(CXX) $^ -shared -o $@

################################################################################
