# -*- makefile -*-
################################################################################
#                                                                               
# Makefile for Walter Dehnen's utils and falcON directories                     
#                                                                               
################################################################################

PWD	:=	$(PWD)/falcON

.PHONY	:	utils falcON
all	:	utils falcON
utils	:
		$(MAKE) -C utils default
ifdef NEMO
falcON	:	utils
		$(MAKE) -C falcON tonemo
else
falcON	:	utils
		$(MAKE) -C falcON default
endif
clean	:
		$(MAKE) -C utils clean
		$(MAKE) -C falcON clean
################################################################################
