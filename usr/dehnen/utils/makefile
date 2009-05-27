# -*- makefile -*-
################################################################################
#
# Makefile for WDutils
#
################################################################################

SHELL	:=	csh

.PHONY		: lib libpgi libpath libsun libicc

default		: lib

# these variables have value 1 if the corresponding compiler is present
# NOTE: we need csh (default for make is sh, which wouldn't work)
__GCC		:= $(shell which gcc    |& wc -w)
__PGI		:= $(shell which pgcc   |& wc -w)
__PTH		:= $(shell which pathcc |& wc -w)
__SUN		:= $(shell which cc     |& wc -w)
__ICC		:= $(shell which icc    |& wc -w)

# GNU compiler (default)
#    used for developing and verifying the code (version 4.3.1)
ifeq ($(__GCC),1)
lib		:
		make -f Makefile COMPILER=gcc
endif

# PGI compiler
#    somewhat tested, seems to work under version 7.2.1
ifeq ($(__PGI),1)
libpgi		:
		make -f Makefile COMPILER=pgi
endif

# Pathscale compiler
#    NOT TESTED
ifeq ($(__PTH),1)
libpath		:
		make -f Makefile COMPILER=pth
endif

# Sun compiler
#    NOT TESTED
ifeq ($(__SUN),1)
libsun		:
		make -f Makefile COMPILER=sun
endif

# Intel compiler
#    NOT TESTED
ifeq ($(__ICC),1)
libicc		:
		make -f Makefile COMPILER=icc
endif

# END