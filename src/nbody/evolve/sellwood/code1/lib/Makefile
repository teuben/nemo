# $NEMO/src/nbody/evolve/sellwood/code1/lib
include $(NEMOLIB)/makedefs

# V1.0	6-jun-97
# V1.1 24-jun-97	added icheck.f measure.f

SRC = bdpres.f fftsna.f grenfn.f  psnset.f vbdass.f \
 blkdat.f        fftsns.f        grget.f         putout.f        vbdgrn.f \
 cenvel.f        findf.f         grnorg.f        slvptf.f        vbdpot.f \
 crash.f         getacc.f        loadup.f        slvptt.f        vbdtab.f \
 dbwrit.f        getgrp.f        markrs.f        step.f          vfthil.f \
 difpot.f        getline.f       masscl.f        stogrp.f        vftsna.f \
 dtrfm1.f        grdset.f        masset.f        stpgrp.f        vftsns.f \
 dtrfm3.f        greenf.f        massgn.f        tmcent.f        weight.f \
 fftcsa.f        greenl.f        newfil.f        tmcgrp.f \
 fftcss.f        greenm.f        opnfil.f        transp.f \
 fftset.f        grenf3.f        poiss.f         unload.f \
 icheck.f  measure.f

INC = ../admin.h ../rjinclude.h


LINKS = admin.h rjinclude.h

OBJ = $(SRC:f=o)

# see somments in ../Makefile on the FFLAGS issue on 64bit
# FFLAGS = -O

LIB = ../galaxy.a

help:
	@echo No help

clean:
	rm -f core *.o *.a 

lib:	$(LINKS) $(OBJ)
	ar ruv $(LIB) $?
	ranlib $(LIB)

links:	$(LINKS)

admin.h: 
	ln -s ../admin.h

rjinclude.h:
	ln -s ../rjinclude.h


$(OBJ): $(INC)

$(SRC): $(INC)
