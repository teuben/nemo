# -*- makefile -*-
################################################################################
#
# makefile for WDutils
#
################################################################################

default			: library
install			: library cleanobj

# -----------------
# Hardware platform
# -----------------
PLATFORM 		:= $(shell uname -m)

# ---------
# MPI stuff
# ---------
MPIDIR			:= $(shell ls | grep mpich)
MPICC			:= $(MPIDIR)/bin/mpicc
MPIINC			:= $(MPIDIR)/include/
IMPI			:= -I$(MPIINC)

# -----------
# directories
# -----------

INC             	:= inc/
SRC             	:= src/
LIB             	:= lib/
LIBT			:= $(LIB).dummy
$(LIBT)			:
			mkdir --parents --mode=755 $(LIB)
			touch $(LIBT)

# -----------------------
# compiler and linker etc
# -----------------------

#DEBUG			:= -DDEBUG

COMPILER		:= gcc#		# GNU C compiler
#COMPILER		:= pgcc#	# problems with linking
#COMPILER		:= icc#		# intel compiler

# -----------------------
# GNU's gcc compiler
ifeq ($(COMPILER),gcc)

CXX			:= $(CPATH)g++
CC			:= $(CPATH)gcc
FC			:= $(CPATH)g77
# CPU dependent optimization flags
# EDIT to match your hardware

ifeq ($(HOSTNAME),virgo)
# on virgo use
CPUFLAGS		:= -march=pentium-m -mpreferred-stack-boundary=4

else

ifeq ($(PLATFORM),i386)
# on i386 use
CPUFLAGS		:= -march=i386 -mpreferred-stack-boundary=4
endif

ifeq ($(PLATFORM),i686)
# on i686 use
CPUFLAGS		:= -march=i686 -mpreferred-stack-boundary=4
endif

ifeq ($(PLATFORM),x86_64)
# on opteron (x86-64) use
CPUFLAGS		:= -march=k8 -mpreferred-stack-boundary=4
endif
endif

ifdef DEBUG
OPTFLAGS		:= -ggdb3 -Wall -Wextra -Wshadow -fPIC $(CPUFLAGS) -rdynamic
else
OPTFLAGS		:= -ggdb3 -Wall -Wextra -Wshadow -O2 -fPIC $(CPUFLAGS) -ffast-math -funroll-loops -fforce-addr -rdynamic
#OPTFLAGS		:= -ggdb3 -O2 -fPIC $(CPUFLAGS) -ffast-math -funroll-loops -fforce-addr -fomit-frame-pointer -rdynamic
endif

else
# -----------------------
# portland group compiler
ifeq ($(COMPILER),pgcc)

CXX			:= pgCC
CC			:= pgcc
FC			:= pgf77

ifdef DEBUG
OPTFLAGS		:= -g -inline_debug_info -fPIC -fastsse --no_using_std --no_exceptions -Msafeptr -Minform=inform -D__PGCC__
else
OPTFLAGS		:= -fast -fastsse -fPIC --no_using_std -Msafeptr -Minform=inform -D__PGCC__
endif


else
# -----------------------
# intel compiler
ifeq ($(COMPILER),icc)

CXX			:= $(CPATH)icpc
CC			:= $(CPATH)icc
FC			:= $(CPATH)ifc
ifdef DEBUG
OPTFLAGS 		:= -g  $(STATIC)
else
#OPTFLAGS		:= -O2 -unroll -tpp7 -xW -i_dynamic -prefetch
OPTFLAGS		:= -O2 -unroll -tpp7 -wd177,1572 -i_dynamic -prefetch  $(STATIC)
endif

else

top			:
			@echo "unknown compiler"
			@echo "please set COMPILER in makefile"

endif
endif
endif

AR			:= ar rc
RL			:= ranlib

# -------------------------
# compiler and linker flags
# -------------------------

CFLAGS			:= $(OPTFLAGS)
CXXFLAGS		:= $(OPTFLAGS)

# -------------------------
# compile and link commands
# -------------------------

MAKE_OBJ		= $(CXX) -c -o $@ $< -I$(INC) $(CXXFLAGS)

# -----------------------
# header dependency lists
# -----------------------

Pi_h			:= $(INC)Pi.h
heap_h			:= $(INC)heap.h
meta_h			:= $(INC)meta.h
inline_h		:= $(INC)inline.h
exception_h		:= $(INC)exception.h
traits_h		:= $(INC)traits.h $(exception_h)
io_h			:= $(INC)io.h $(exception_h)
tupel_h			:= $(INC)tupel.h $(INC)tupel.cc $(traits_h)
memory_h		:= $(INC)memory.h $(exception_h) $(traits_h) $(inline_h)
numerics_h		:= $(INC)numerics.h $(inline_h) $(memory_h) $(tupel_h)
spline_h		:= $(INC)spline.h $(numerics_h) $(memory_h)
random_h		:= $(INC)random.h $(inline_h) $(Pi_h) $(traits_h)
WDMath_h		:= $(INC)WDMath.h $(Pi_h) $(exception_h) $(traits_h) \
				$(inline_h) 
parallel_h		:= $(INC)parallel.h $(traits_h) $(memory_h) $(tupel_h)

# -----------------------
# source dependency lists
# -----------------------

meta_cc			:= $(SRC)meta.cc $(meta_h)
exception_cc		:= $(SRC)exception.cc $(exception_h)
io_cc			:= $(SRC)io.cc $(io_h)
numerics_cc		:= $(SRC)numerics.cc $(numerics_h) $(WDMath_h)
random_cc		:= $(SRC)random.cc $(random_h) $(numerics_h) \
				$(exception_h)
WDMath_cc		:= $(SRC)WDMath.cc $(WDMath_h) $(inline_h)
parallel_cc		:= $(SRC)parallel.cc $(parallel_h) $(meta_h) $(io_h)

# ---------------
# library modules
# ---------------

$(LIB)exception.o:	$(exception_cc) $(LIBT) makefile
			$(MAKE_OBJ)
$(LIB)io.o:		$(io_cc) $(LIBT) makefile
			$(MAKE_OBJ)
$(LIB)meta.o:		$(meta_cc) $(LIBT) makefile
			$(MAKE_OBJ)
$(LIB)numerics.o:	$(numerics_cc) $(LIBT) makefile
			$(MAKE_OBJ)
$(LIB)random.o:		$(random_cc) $(LIBT) makefile
			$(MAKE_OBJ)
$(LIB)WDMath.o:		$(WDMath_cc) $(LIBT) makefile
			$(MAKE_OBJ)

WDutils_objs	:=	$(LIB)exception.o $(LIB)io.o $(LIB)meta.o \
			$(LIB)numerics.o $(LIB)random.o $(LIB)WDMath.o

# -------
# library
# -------

$(LIB)libWDutils.a:	$(WDutils_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)libWDutils.so:	$(WDutils_objs)
			$(CXX) $^ -shared -o $@ $(CPUFLAGS)

library		:	$(LIB)libWDutils.a $(LIB)libWDutils.so

ifeq ($(MPIDIR),mpich)
# -------------------
# MPI library modules
# -------------------

$(LIB)parallel.o:	$(parallel_cc) $(makefiles)
			$(MAKE_OBJ) $(IMPI)

WDutilsP_objs	:=	$(LIB)exception.o $(LIB)io.o $(LIB)meta.o \
			$(LIB)numerics.o $(LIB)random.o $(LIB)WDMath.o \
			$(LIB)parallel.o

# -----------
# MPI library
# -----------

$(LIB)libWDutilsMPI.a:	$(WDutilsP_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)libWDutilsMPI.so:	$(WDutilsP_objs)
			$(CXX) $^ -shared -o $@ $(CPUFLAGS)

library		:	$(LIB)libWDutils.a $(LIB)libWDutils.so \
			$(LIB)libWDutilsMPI.a $(LIB)libWDutilsMPI.so

else

library		:	$(LIB)libWDutils.a $(LIB)libWDutils.so

endif

# --------
# cleaning
# --------

.PHONY		: 	cleanbackup cleandirs clean cleanall
cleanbackup	:
			rm -f *~ $(INC)*~ $(SRC)*~
cleandirs	:
			rm -rf $(LIB)

clean		:	cleandirs

cleanall	:	cleandirs cleanbackup
