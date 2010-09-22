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

# ----------------
# development flag
# ----------------

ifndef Standard

ifeq ($(USER),wd11)
WDutilsDevel	:= -DWDutilsDevel
endif

ifeq ($(USER),walter)
WDutilsDevel	:= -DWDutilsDevel
endif

endif

# ---------
# libraries
# ---------

.PHONY		: libgcc libpgi libpath libsun libicc test

ifdef WDutilsDevel
default		: libgcc testsuite
else
default		: libgcc
endif

# these variables have value "1" or ("       1" on a MAC) if the
# corresponding compiler is present
# NOTE: we need csh (default for make is sh, which wouldn't work)
__GCC		:= $(shell which gcc    |& wc -w)
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
		make -f makeall COMPILER=gcc WDutilsDevel=$(WDutilsDevel) library
ifdef WDutilsDevel
testsuite	:
		make -f makeall COMPILER=gcc WDutilsDevel=$(WDutilsDevel) testsuite
endif
endif

# PGI compiler
#    somewhat tested, seems to work under version 7.2.1
ifeq ($(__PGI),$(VALID))
libpgi		:
		make -f makeall COMPILER=pgi WDutilsDevel=$(WDutilsDevel) library
endif

# Pathscale compiler
#    NOT TESTED
ifeq ($(__PTH),$(VALID))
libpath		:
		make -f makeall COMPILER=pth WDutilsDevel=$(WDutilsDevel) library
endif

# Sun compiler
#    NOT TESTED
ifeq ($(__SUN),$(VALID))
libsun		:
		make -f makeall COMPILER=sun WDutilsDevel=$(WDutilsDevel) library
endif

# Intel compiler
#    NOT TESTED
ifeq ($(__ICC),$(VALID))
libicc		:
		make -f makeall COMPILER=icc WDutilsDevel=$(WDutilsDevel) library
ifdef WDutilsDevel
testsuiteicc	:
		make -f makeall COMPILER=icc WDutilsDevel=$(WDutilsDevel) testsuite
endif
endif

# ---------------------
# doxygen documentation
# ---------------------

doxu		:
		make -f makeall COMPILER=gcc WDutilsDevel=$(WDutilsDevel) doxu

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