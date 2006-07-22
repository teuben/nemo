# -*- makefile -*-
################################################################################
#
# makefile for WDutils
#
################################################################################

default			: libs
install			: libs cleanobj

# ------------------
# Hardware platform
# ------------------
PLATFORM := $(shell uname -m)

# -----------
# directories
# -----------

INC             	:= inc/
SRC             	:= src/
LIB             	:= lib/$(PLATFORM)/
.PHONY          	: dirs
dirs            	:
			./make.dir lib    755
			./make.dir $(LIB) 755

# -----------------------
# compiler and linker etc
# -----------------------

COMPILER		:= gcc
#COMPILER		:= pgcc#         problems with linking
#COMPILER		:= icc


ifeq ($(COMPILER),gcc)
# GNU's gcc compiler

CXX			:= $(CPATH)g++
CC			:= $(CPATH)gcc
FC			:= $(CPATH)g77
# CPU dependent optimization flags
# EDIT to match your hardware

ifeq ($(PLATFORM),i386)
# on i386 use
CPUFLAGS		:= -march=i386 -mpreferred-stack-boundary=4
endif

ifeq ($(PLATFORM),i686)
# on i686 use
CPUFLAGS		:= -march=athlon-xp -mpreferred-stack-boundary=4
endif

ifeq ($(PLATFORM),x86_64)
# on opteron (x86-64) use
CPUFLAGS		:= -mcpu=k8 -mpreferred-stack-boundary=4
endif

OPTFLAGS		:= -ggdb -O2 -fPIC $(CPUFLAGS) -funroll-loops -fforce-addr -fomit-frame-pointer -rdynamic

else
ifeq ($(COMPILER),pgcc)
# portland groups compiler

CXX			:= pgCC
CC			:= pgcc
FC			:= pgf77
#OPTFLAGS		:= -fast -Minline=levels:10 -Mcache_align --no_exceptions -Minform=inform -D__PGCC__
ifdef DEBUG
OPTFLAGS		:= -g -inline_debug_info -O2 -fastsse --no_using_std --no_exceptions -Msafeptr -Minform=inform -D__PGCC__
else
OPTFLAGS		:= -O2 -fastsse --no_using_std --no_exceptions -Msafeptr -Minform=inform
-D__PGCC__
endif


else
ifeq ($(COMPILER),icc)
# intel compiler

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
inline_h		:= $(INC)inline.h
inline_io_h		:= $(INC)inline_io.h
exception_h		:= $(INC)exception.h
traits_h		:= $(INC)traits.h $(exception_h)
tupel_h			:= $(INC)tupel.h $(INC)tupel.cc $(traits_h)
memory_h		:= $(INC)memory.h $(exception_h) $(traits_h) $(inline_h)
numerics_h		:= $(INC)numerics.h $(inline_h) $(memory_h) $(tupel_h)
spline_h		:= $(INC)spline.h $(numerics_h) $(memory_h)
random_h		:= $(INC)random.h $(inline_h) $(Pi_h) $(traits_h)
WDMath_h		:= $(INC)WDMath.h $(Pi_h) $(exception_h) $(traits_h) \
				$(inline_h) 

# -----------------------
# source dependency lists
# -----------------------

exception_cc		:= $(SRC)exception.cc $(exception_h)
numerics_cc		:= $(SRC)numerics.cc $(numerics_h) $(WDMath_h)
random_cc		:= $(SRC)random.cc $(random_h) $(numerics_h) \
				$(exception_h)
WDMath_cc			:= $(SRC)WDMath.cc $(WDMath_h) $(inline_h)

# ---------------
# library modules
# ---------------

$(LIB)exception.o:	$(exception_cc) makefile
			./make.dir lib    755
			./make.dir $(LIB) 755
			$(MAKE_OBJ)
$(LIB)numerics.o:	$(numerics_cc) makefile
			$(MAKE_OBJ)
$(LIB)random.o:		$(random_cc) makefile
			$(MAKE_OBJ)
$(LIB)WDMath.o:		$(WDMath_cc) makefile
			$(MAKE_OBJ)

WDutils_objs	:=	$(LIB)exception.o $(LIB)numerics.o $(LIB)random.o \
			$(LIB)WDMath.o

# -------
# library
# -------

$(LIB)libWDutils.a:	$(WDutils_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)libWDutils.so:	$(WDutils_objs)
			$(CXX) $^ -shared -o $@ $(CPUFLAGS)

libs			: $(LIB)libWDutils.a $(LIB)libWDutils.so

# --------
# cleaning
# --------

.PHONY		: 	cleanbackup cleanlib cleanobj clean
cleanbackup:
			rm -f *~ $(INC)*~ $(INC)*~ $(SRC)*~
cleanobj:
			rm -f $(LIB)*.o
cleanlib:
			rm -f $(LIB)*
clean:			cleanbackup cleanlib
