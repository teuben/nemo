# -*- makefile -*-
################################################################################
#
# Makefile for WDutils
#
################################################################################

SHELL		:= csh

# -----------
# directories
# -----------

INC             := inc/
SRC             := src/

# --------------------------
# gcc compiler and testsuite
# --------------------------
GCC		:= gcc
GPP		:= g++

ifeq ($(VENDOR),apple)
ifeq ($(USER),walter)
WDutilsTEST	:= -DWDutilsTEST
# on my MAC, I'v installed gcc 4.4.1 under /sw.
# it's called gcc-4 and g++-4
GCC		:= gcc-4
GPP		:= g++-4
endif
else
ifeq ($(USER),wd11)
WDutilsTEST	:= -DWDutilsTEST
endif
endif

# ---------
# libraries
# ---------

.PHONY		: libgcc libpgi libpath libsun libicc test

default		: libgcc

# these variables have value "1" or ("       1" on a MAC) if the
# corresponding compiler is present
# NOTE: we need csh (default for make is sh, which wouldn't work)
__GCC		:= $(shell which $(GCC) |& wc -w)
__PGI		:= $(shell which pgcc   |& wc -w)
__PTH		:= $(shell which pathcc |& wc -w)
__SUN		:= $(shell which cc     |& wc -w)
__ICC		:= $(shell which icc    |& wc -w)
# this is the expected value for available compilers
VALID		:= $(shell which make   |& wc -w)


# GNU compiler (default)
#    used for developing and verifying the code (version 4.3 or higher)
ifeq ($(__GCC),$(VALID))
libgcc		:
		make -f makeall COMPILER=gcc GCC=$(GCC) GPP=$(GPP) library
ifdef WDutilsTEST
test		:
		make -f makeall COMPILER=gcc GCC=$(GCC) GPP=$(GPP) test
endif
endif

# PGI compiler
#    somewhat tested, seems to work under version 7.2.1
ifeq ($(__PGI),$(VALID))
libpgi		:
		make -f makeall COMPILER=pgi library
endif

# Pathscale compiler
#    NOT TESTED
ifeq ($(__PTH),$(VALID))
libpath		:
		make -f makeall COMPILER=pth library
endif

# Sun compiler
#    NOT TESTED
ifeq ($(__SUN),$(VALID))
libsun		:
		make -f makeall COMPILER=sun library
endif

# Intel compiler
#    NOT TESTED
ifeq ($(__ICC),$(VALID))
libicc		:
		make -f makeall COMPILER=icc library
endif

# ---------------------
# doxygen documentation
# ---------------------

doxu		:
		make -f makeall COMPILER=gcc doxu

# --------
# cleaning
# --------

.PHONY		: 	cleanbackup cleandirs clean cleanall
cleanbackup	:
			rm -f *~ $(INC)*~ $(SRC)*~
cleandirs	:
			rm -rf lib lib.pgi lib.path lib.sun lib.icc

clean		:	cleandirs

cleanall	:	cleandirs cleanbackup

# END