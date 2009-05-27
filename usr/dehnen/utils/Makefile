# -*- makefile -*-
################################################################################
#
# Makefile for WDutils
#
################################################################################

default			: library
install			: library cleanobj

# ----------------
# Compiler related
# ----------------

# GNU compiler
ifeq ($(COMPILER),gcc)
-include make.gcc
endif

# Portland group compiler
ifeq ($(COMPILER),pgi)
-include make.pgi
endif

# Pathscale compiler
ifeq ($(COMPILER),pth)
-include make.pth
endif

# Sun compiler
ifeq ($(COMPILER),sun)
-include make.sun
endif

# Intel compiler
ifeq ($(COMPILER),icc)
-include make.icc
endif


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
LIBT			:= $(LIB).dummy
$(LIBT)			:
			mkdir --parents --mode=755 $(LIB)
			touch $(LIBT)

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

AR			:= ar rc
RL			:= ranlib

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
