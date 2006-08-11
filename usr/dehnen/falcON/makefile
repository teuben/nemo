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

ifdef DPROPER
LIBNAME			:= falcON.P
else
LIBNAME			:= falcON
endif

#
# 2 set compiler
# this should be set by a configure file
#
COMPILER		:= gcc
#COMPILER		:= pgcc#	 problems with linking
#COMPILER		:= icc

#
# 3 activate options (some only apply to proprietary and private versions)
#
ifdef DPROPER
ifeq ($(COMPILER),gcc)
# 4.1 use SSE instructions in force computation (proprietary only)
DSSE 			:= -DfalcON_SSE
endif
endif

ifdef DPROPER
# 3.2 full sph support (development version)
DSPH			:= -DfalcON_SPH
endif
ifdef DWALTER
# 3.3 MPI parallel code (under development)
#DMPI			:= -DfalcON_MPI
endif

# 3.4 individual softening lengths (default: non-adaptive)
DSOFT 			:= -DfalcON_INDI# -DfalcON_ADAP

# 3.5 precision of floating point numbers (default: 32 bit)
DPRECISION		:= -DfalcON_SINGLE#       32 bit    default
#DPRECISION		:= -DfalcON_DOUBLE#       64 bit


# 3.6 compile with little optimisation
#DEBUG			:= -DEBUG -DDEBUG

# 3.7 use static library
#     NOTE: this will break make since executables are no re-made if library is
ifdef DEBUG
#STATIC		:= -static
endif

#
# 4 symbolic links
#

ifdef DPROPER
ifdef DWALTER
ifdef DMPI
links			:= links_pro links_sph links_walter links_pjm links_mpi
else
links			:= links_pro links_sph links_walter links_pjm
endif
else
links			:= links_pro links_sph
endif
endif

#
# 5 global pseudo targets
#
.PHONY			: manip acc man utils

ifdef NEMO
default			: exe manip acc doxu
ifdef DPROPER
# proprietary version: don't copy executables and manipulators
# NOTE: this means we have to set the PATH and MANIPPATH appropriately
tonemo			: exe manip acc doxu
				cp $(LIB)/*.so $(NEMOLIB)
			      - cp $(ACC)/* $(NEMOOBJ)/acc/
else
# public version: copy executables, manipulators, and accelerations
tonemo			: exe manip acc doxu
				cp $(LIB)/*.so $(NEMOLIB)
			      - cp $(BIN)/* $(NEMOBIN)
			      - cp $(MAN)/* $(NEMOOBJ)/manip/
			      - cp $(ACC)/* $(NEMOOBJ)/acc/
endif
install			: tonemo man
else
default			: exe
install			: default
endif

#
# 6 include specific makefiles
#


# 6.1 makedefs: defines compiler options etc
-include makedefs
# 6.2 makepub: all targets for the public version
-include makepub
ifdef DPROPER
# 6.3 proper/make: additional targets for proprietary version
-include proper/make
endif

ifdef DSPH
# 6.4 sph/make: additional targets for SPH supporting version
-include sph/make
endif

ifdef DMPI
# 6.5 parallel/make: additional targets for MPI parallel version
-include parallel/make
endif

ifdef DWALTER
# 6.6 additional targets for any private sandboxes
-include walter/make
-include pjm/make

#
# 7 pseudo targets to make various things ...
#

# 7.1 user guides
guide			:
			$(MAKE) -C tex default
endif

# 7.1 utils library
utils			:
			ln -sf ../utils
			(cd inc; ln -sf ../utils/inc utils)
			$(MAKE) -C utils default

# 7.2 falcON library
# NOTE: $(falcON) defined in makedefs, depending on value of STATIC
falcON_objs	=	$(public_objs) $(proper_objs) $(sph_objs) $(parallel_objs)

$(LIB)lib$(LIBNAME).a:	$(falcON_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)lib$(LIBNAME).so:	$(falcON_objs)
			$(CXX) $^ -shared -o $@

library		: 	utils $(LIBT) $(falcON) $(BFCT)

# 7.3 falcON executables
exe		:	library $(exe_pub) $(exe_pro) $(exe_sph) $(exe_walter) $(exe_pjm)

ifdef NEMO

# 7.4 falcON manipulators
manip		:	utils $(MANT) $(manip_pub) $(manip_pro) $(manip_walter)

# 7.5 falcON accelerations
acc		:	$(ACCT) $(acc_pub) $(acc_pro) $(acc_walter) $(acc_pjm)	

# 7.6 copy falcON man pages to nemo
ifdef DPROPER
man		:
else
man		:
			- cp -p man/man1/*.? $(NEMO)/man/man1
			- cp -p man/man3/*.? $(NEMO)/man/man3
			- cp -p man/man5/*.? $(NEMO)/man/man5
endif
endif

# 7.7 several clean options
# 7.7.1 clean editor backup files
cleanbackup	:
			rm -f *~ */*~ */*/*~  */*/*/*~
# 7.7.2 clean objects (if no longer needed, since either in library or executable)
cleanobjs	:
			rm -rf $(LIB)*.o

# 7.7.3 cleans everything made in releases
cleanrelease	:
			rm -rf $(BIN) $(MAN) $(ACC)
			rm -rf dox/html

# 7.7.4 clean == cleanrelease
clean		:	cleanrelease

# 7.7.5 cleanall == cleandirs & cleanbackup and remove user guides
cleanall	:	cleanrelease cleanbackup
			rm doc/user*.pdf

# 7.7.6 everything to be done before we can checkin a public version
beforecheckin	:	cleanbackup
			rm doc/user_guide_proper.pdf
			$(MAKE) -C tex public

################################################################################
