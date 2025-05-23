# -*- makefile -*-
################################################################################
#                                                                               
# Makefile for the falcON project                                               
#                                                                               
# this makefile should work for the public, proprietary, and private versions   
#                                                                               
################################################################################

#
# 1 name of library
#

ifndef LIBNAME
ifdef DPROPER
LIBNAME			:= falcON.P
else
LIBNAME			:= falcON
endif
endif

#
# 2 set compiler
# this should be set by a configure file
#

ifeq ($(COMPILER),)
COMPILER		:= gcc#		# GNU compiler g++
endif
#COMPILER		:= icc#		# intel compiler icpc
#COMPILER		:= clang#	# llvm compiler clang++


#
# 3 activate options (some only apply to proprietary and private versions)
#
ifdef DPROPER
ifeq ($(COMPILER),gcc)
# 3.1 use SSE instructions in force computation (proprietary only)
#
# NOTE with gcc 4.1.2 this does not work anymore and has to be commented out
#
#DSSE 			:= -DfalcON_SSE
endif
endif

# 3.3 individual softening lengths (default: non-adaptive)
DSOFT 			:= # -DfalcON_ADAP

# 3.4 precision of floating point numbers (default: 32 bit)
DPRECISION		:= -DfalcON_SINGLE#       32 bit    default
#DPRECISION		:= -DfalcON_DOUBLE#       64 bit


# 3.5 compile with little optimisation
#DEBUG			:= -DEBUG -DDEBUG

# 3.6 use static library
#     NOTE: this will break make since executables are no re-made if library is
ifdef DEBUG
#STATIC		:= -static
endif

#
# 4 global pseudo targets
#
.PHONY			: manip acc man

ifdef NEMO
default			: links exe manip acc doxu
ifdef DPROPER
# proprietary version: don't copy executables and manipulators
# NOTE: this means we have to set the PATH and MANIPPATH appropriately
tonemo			: acc doxu man
				cp $(ACC)/* $(NEMOOBJ)/acc/
else
# public version: copy executables, manipulators, and accelerations
tonemo			: exe manip acc doxu man
				cp $(BIN)/* $(NEMOBIN)
			      - cp $(MAN)/* $(NEMOOBJ)/manip/
			      - cp $(ACC)/* $(NEMOOBJ)/acc/
endif
install			: tonemo
else
default			: exe
install			: default
endif

#
# 5 include specific makefiles
#


# 5.1 makedefs: defines compiler options etc
-include makedefs

# 5.2 makepub: all targets for the public version
-include makepub

ifndef Standard

# 5.3 proper/make: additional targets for proprietary version
ifdef DPROPER
-include proper/make
endif

# 5.4 sph/make: additional targets for SPH supporting version
ifdef DSPH
-include sph/make
endif

# 5.5 additional targets for any private sandboxes
ifdef DWALTER
-include walter/make
endif

endif

#
# 6 pseudo targets to make various things ...
#

# 6.1 user guides
ifdef DWALTER
guide			:
			$(MAKE) -C tex default
endif

# 6.1 symbolic links
utils			:
			ln -sf ../utils
inc/utils		: utils
			(cd inc; ln -sf ../utils/inc utils)
links			:=	$(links_pro) utils inc/utils $(links_sph) $(links_walter) $(links_pjm) $(links_mpi)

links			:	$(links)

# 6.1 utils library
utils/lib/libWDutils.a	: utils
			$(MAKE) -C utils
libutils		: utils utils/lib/libWDutils.a

# 6.3 falcON library
# NOTE: $(falcON) defined in makedefs, depending on value of STATIC
falcON_objs	=	$(public_objs) $(proper_objs) $(sph_objs)

$(LIB)lib$(LIBNAME).a:	$(falcON_objs)
			$(AR) $@ $?
			$(RL) $@

ifeq ($(NO_SHARED_LIBS),0)
$(LIB)lib$(LIBNAME).so:	$(falcON_objs)
			$(CXX) $^ -L$(NEMO)/usr/dehnen/utils/lib -lWDutils $(LNEMO) -shared -o $@
ifeq ($(OS),Darwin)
			install_name_tool -id $(NEMO)/usr/dehnen/falcON/lib/lib$(LIBNAME).so $(LIB)lib$(LIBNAME).so
endif
else
$(LIB)lib$(LIBNAME).so:
			@echo Skipping lib$(LIBNAME).so:
endif

library		: 	links $(LIBT) $(falcON) $(BFCT)

# 6.4 falcON executables
exe		:	libutils library $(exe_pub) $(exe_pro) $(exe_sph) $(exe_walter) $(exe_pjm) $(exe_mpi)

ifdef NEMO

# 6.5 falcON manipulators
manip		:	links $(MANT) $(manip_pub) $(manip_pro) $(manip_sph) $(manip_mpi) $(manip_walter)

# 6.6 falcON accelerations
acc		:	$(ACCT) $(acc_pub) $(acc_pro) $(acc_sph) $(acc_walter) $(acc_pjm)	

# 6.7 copy falcON man pages to nemo
ifdef DPROPER
man		:
else
man	        :
			@echo "no man page sync to NEMO, assuming NEMO has the latest"
			@echo "use 'make man_orig' to force, but check git logs"
man_orig		:
			- cp -p man/man1/*.? $(NEMO)/man/man1
			- cp -p man/man3/*.? $(NEMO)/man/man3
			- cp -p man/man5/*.? $(NEMO)/man/man5
endif
endif

# 6.8 several clean options
# 6.8.1 clean editor backup files
cleanbackup	:
			rm -f *~ */*~ */*/*~  */*/*/*~
# 6.8.2 clean objects (if no longer needed, since either in library or executable)
cleanobjs	:
			rm -rf $(OBJ)*.o $(LIB)*.so $(LIB)*.a

# 6.8.3 cleans symbolic links
cleanlinks	:	
			rm $(links)

# 6.8.3 cleans everything made in releases
cleanrelease	:
			rm -rf $(BIN) $(LIB) $(MAN) $(ACC)
			rm -rf dox/html

# 6.8.4 clean == cleanrelease
clean		:	cleanrelease

# 6.8.5 cleanall == cleandirs & cleanbackup and remove user guides
cleanall	:	cleanrelease cleanbackup
			rm doc/user*.pdf

# 6.8.6 everything to be done before we can checkin a public version
beforecheckin	:	cleanbackup
			$(MAKE) -C tex public

################################################################################
