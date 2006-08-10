# -*- makefile -*-
################################################################################
#                                                                               
# Makefile for the falcON project                                               
#                                                                               
################################################################################
.PHONY			: manip acc man utils

ifdef NEMO
default			: exe manip acc doxu
tonemo			: exe manip acc doxu
				cp $(LIB)/libfalcON.so $(NEMOLIB)
			      - cp $(BIN)/* $(NEMOBIN)
			      - cp $(MAN)/* $(NEMOOBJ)/manip/
			      - cp $(ACC)/* $(NEMOOBJ)/acc/
install			: tonemo man
else
default			: exe
install			: default
endif


#
# this should be set by a configure file
#
ifeq ($(USER),wd11)
ifeq ($(PWD),/home/wd11/Code/falcON)
DWALTER			:= -DWalter
endif
endif

ifdef DWALTER
DPROPER 		:= -DfalcON_PROPER
DSPH			:= -DfalcON_SPH#    comment out to disable SPH code     
##DMPI			:= -DfalcON_MPI#    comment out to disable MPI code     
else
#DPROPER 		:= -DfalcON_PROPER# comment out to obtain public version
endif

#
# this should be set by a configure file
#
#DPRECISION		:= -DfalcON_SINGLE#       32 bit    default
#DPRECISION		:= -DfalcON_DOUBLE#       64 bit

#
# this should be set by a configure file
#
COMPILER		:= gcc
#COMPILER		:= pgcc#	 problems with linking
#COMPILER		:= icc

DSOFT 			:= -DfalcON_INDI# -DfalcON_ADAP
                                        #   individual softening lengths        

ifdef DPROPER
ifeq ($(COMPILER),gcc)
DSSE 			:= -DfalcON_SSE#    comment out to disable SSE code     
endif
endif

##DEBUG			:= -DEBUG -DDEBUG

ifdef DEBUG
##STATIC		:= -static
endif

include makedefs
include makepub
ifdef DPROPER
include proper/make
endif

ifdef DSPH
include sph/make
endif

ifdef DMPI
include parallel/make
endif

ifdef DWALTER
#-include walter/make
#-include pjm/make
endif

#
# list of pseudo targets to make various things ...
#

# user guides
guide			:
			$(MAKE) -C tex default

# utils
utils			:
			ln -sf ../utils
			(cd inc; ln -sf ../utils/inc utils)
			$(MAKE) -C utils default

# falcON library
falcON_objs	=	$(public_objs) $(proper_objs) $(sph_objs) $(parallel_objs)

$(LIB)libfalcON.a:	$(falcON_objs)
			$(AR) $@ $?
			$(RL) $@
$(LIB)libfalcON.so:	$(falcON_objs)
			$(CXX) $^ -shared -o $@

# NOTE: $(falcON) defined in makedefs, depending on value of STATIC
library		: 	utils $(LIBT) $(falcON)

# falcON executables
exe		:	library $(exe_pub) $(exe_pro) $(exe_sph) $(exe_walter)

ifdef NEMO

# falcON manipulators
manip		:	utils $(MANT) $(manip_pub) $(manip_pro) $(manip_walter)

# falcON accelerations
acc		:	$(ACCT) $(acc_pub) $(acc_pro)	

# copy falcON man pages to nemo
man		:
			- cp -p man/man1/*.? $(NEMO)/man/man1
			- cp -p man/man3/*.? $(NEMO)/man/man3
			- cp -p man/man5/*.? $(NEMO)/man/man5
endif

# clean editor backup files
cleanbackup	:
			rm -f *~ $(INC)*~ $(INC)*/*~ $(SRC)*/*~ $(SRC)*/*/*~
# clean objects (if no longer needed, since either in library or executable)
cleanobjs	:
			rm -rf $(LIB)*.o

# cleans everything made in releases
cleanrelease	:
			rm -rf $(BIN) $(MAN) $(ACC)
			rm -rf dox/html

# clean == cleanrelease
clean		:	cleanrelease

# cleanall == cleandirs & cleanbackup and remove user guides
cleanall	:	cleanrelease cleanbackup
			rm doc/user*.pdf

# everything to be done before we can checkin a public version
beforecheckin	:	cleanbackup
			rm doc/user_guide_proper.pdf
			$(MAKE) -C tex public

################################################################################
