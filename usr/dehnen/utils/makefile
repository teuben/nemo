# -*- makefile -*-
################################################################################
#
# Makefile for WDutils
#
################################################################################

SHELL		:= csh

default		: library

# ----------------
# Compiler related
# ----------------

ifeq ($(FALCONEXT),.gcc)
  COMPILER := gcc
endif
ifeq ($(FALCONEXT),.clang)
  COMPILER := clang
endif
ifeq ($(FALCONEXT),.icc)
  COMPILER := icc
endif
ifeq ($(COMPILER),)
  COMPILER := gcc
endif

# NOTE: EXT is re-defined in the make.$(EXT) files to include a leading '.'

include make.$(COMPILER)

# GNU compiler
ifeq ($(COMPILER),gcc)
#SAFE_DIRECT	:= -fno-associative-math
#CGAL_FLAGS	:= -frounding-math
endif

# -----------
# directories
# -----------

MODE			:= 755
PUBLICMODE		:= 777

INC             	:= inc/
SRC             	:= src/
OBJ			:= obj/
OBJT			:= $(OBJ).dummy
$(OBJT)			:
			mkdir -p -m 755 $(OBJ)
			touch $(OBJT)

LIB			:= lib/
LIBT			:= $(LIB).dummy
$(LIBT)			:
			mkdir -p -m 755 $(LIB)
			touch $(LIBT)

# ---------
# libraries
# ---------

INCLUDE			:= -I$(INC)
LINKLIBS		+= $(CXXLIBS)

# math lib
ifndef MATHLIB
MATHLIB			:= -lm
endif

# WDutils lib
INCLUDE			:= -I$(INC)
LUTIL			:= -L$(LIB) -lWDutils

# -------------------------
# compile and link commands
# -------------------------

LINKLIBS		+= $(MATHLIB)
LIBRARIES		:= $(LUTIL) $(LINKLIBS)
MAKE			= $(CXX) -o $@ $< $(INCLUDE) $(CXXFLAGS)
MAKE_OBJ		= $(MAKE) -c
MAKE_EXE		= $(MAKE) $(LIBRARIES)

# -----------------------
# header dependency lists
# -----------------------

$(INC)cachesize.h	: make_cachesize
			tcsh ./make_cachesize > $@
cachesize_h		:= $(INC)cachesize.h
Pi_h			:= $(INC)Pi.h
timer_h			:= $(INC)timer.h
inline_h		:= $(INC)inline.h
exception_h		:= $(INC)exception.h
meta_h			:= $(INC)meta.h $(exception_h)
traits_h		:= $(INC)traits.h $(exception_h)
memory_h		:= $(INC)memory.h $(cachesize_h) $(traits_h) $(inline_h)
io_h			:= $(INC)io.h $(exception_h) $(traits_h)
vector_h		:= $(INC)vector.h $(INC)tupel.h $(INC)tupel.cc
radix_h			:= $(INC)radix.h $(traits_h) $(memory_h)
sse_h			:= $(INC)sse.h $(cachesize_h) $(exception_h) $(meta_h) \
				$(memory_h)
heap_h			:= $(INC)heap.h $(memory_h)
numerics_h		:= $(INC)numerics.h $(inline_h) $(memory_h)
spline_h		:= $(INC)spline.h $(numerics_h) $(memory_h)
random_h		:= $(INC)random.h $(inline_h) $(Pi_h) $(traits_h)
WDMath_h		:= $(INC)WDMath.h $(Pi_h) $(exception_h) $(traits_h) \
				$(inline_h) $(meta_h)
geometry_h		:= $(INC)geometry.h $(sse_h) $(vector_h) \
				$(INC)geometry.tcc
periodic_h		:= $(INC)periodic.h $(geometry_h) $(inline_h)
matr33_h		:= $(INC)matr33.h $(exception_h) $(inline_h)
parse_h			:= $(INC)parse.h

# -----------------------
# source dependency lists
# -----------------------

exception_cc		:= $(SRC)exception.cc $(exception_h)
io_cc			:= $(SRC)io.cc $(io_h) $(memory_h)
numerics_cc		:= $(SRC)numerics.cc $(numerics_h) $(WDMath_h)
random_cc		:= $(SRC)random.cc $(random_h) $(numerics_h) \
				$(exception_h)
WDMath_cc		:= $(SRC)WDMath.cc $(WDMath_h) $(inline_h)

# ---------------
# library modules
# ---------------

$(OBJ)exception.o:	$(exception_cc) $(OBJT) $(makefiles)
			$(MAKE_OBJ)
LIBOBJECTS	+=	$(OBJ)exception.o
$(OBJ)io.o:		$(io_cc) $(OBJT) $(makefiles)
			$(MAKE_OBJ)
LIBOBJECTS	+=	$(OBJ)io.o
$(OBJ)numerics.o:	$(numerics_cc) $(OBJT) $(makefiles)
			$(MAKE_OBJ)
LIBOBJECTS	+=	$(OBJ)numerics.o
$(OBJ)random.o:		$(random_cc) $(OBJT) $(makefiles)
			$(MAKE_OBJ)
LIBOBJECTS	+=	$(OBJ)random.o
$(OBJ)WDMath.o:		$(WDMath_cc) $(OBJT) $(makefiles)
			$(MAKE_OBJ)
LIBOBJECTS	+=	$(OBJ)WDMath.o

# -------
# library
# -------

AR			:= ar rc
RL			:= ranlib

$(LIB)libWDutils.a	:	$(LIBOBJECTS)
				$(AR) $@ $?
				$(RL) $@
$(LIB)libWDutils.so	:	$(LIBOBJECTS)
				$(CXX) $^ $(LDFLAGS) $(LINKLIBS) -o $@

library			:	$(LIBT) $(LIB)libWDutils.a \
					$(LIB)libWDutils.so

# -------------
# test routines
# -------------

test/TestSSE	:	test/TestSSE.cc $(sse_h) $(memory_h) $(makefiles)
			$(MAKE_EXE)
TESTSUITE	+=	test/TestSSE

test/TestGeometry:	test/TestGeometry.cc $(geometry_h) $(makefiles)
			$(MAKE_EXE)
TESTSUITE	+=	test/TestGeometry

test/TestMatrix	:	test/TestMatrix.cc $(matr33_h) $(makefiles)
			$(MAKE_EXE)
TESTSUITE	+=	test/TestMatrix

test/TestLoad	:	test/TestLoad.cc $(load_h) $(exception_h) $(OBJ)load.o
			$(MAKE_EXE) $(OBJ)load.o
TESTSUITE	+=	test/TestLoad

testsuite	:	library $(TESTSUITE)


# ---------------------
# doxygen documentation
# ---------------------

dox/html/index.html:	$(Pi_h) $(heap_h) $(meta_h) $(inline_h) $(exception_h) \
			$(sse_h) $(traits_h) $(io_h) $(vector_h) $(memory_h) \
			$(numerics_h) $(spline_h) $(radix_h) $(random_h) \
			$(WDMath_h)
			doxygen utils.dox

doxu		:	dox/html/index.html

# ------------------
# collective targets
# ------------------

.PHONY		: 	utils cleanbackup cleandirs clean cleanall
testcompiler	:	
			$(CXX) -v
cleanbackup	:
			rm -f *~ $(INC)*~ $(SRC)*~ test/*~
cleandirs	:
			rm -rf obj  lib
clean		:	cleandirs
cleanall	:	cleandirs cleanbackup

# END
