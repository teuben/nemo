# $NEMO/src/nbody/evolve/sellwood/code1
include $(NEMOLIB)/makedefs
#

INC = admin.h  rjinclude.h
SRC = galaxy.f

help:
	@echo GALAXY program, written by Jerry Sellwood

#   on 64bit machines, -O will coredump, -O2 is ok, adding -g also
#   we'll inherit FFLAGS from the makedefs file instead
#FFLAGS = -O2

all:	galaxy.a galaxy

clean:
	rm -f $(LIB) galaxy *.o *.a core galaxy.res galaxy.fin galaxy.lis
	(cd lib; make clean)


galaxy: galaxy.a $(SRC) $(INC)
	$(FC) $(FFLAGS) -o galaxy galaxy.f galaxy.a


galaxy.a: $(INC)
	(cd lib; make lib)

install:	galaxy
	cp $? $(NEMOBIN)
