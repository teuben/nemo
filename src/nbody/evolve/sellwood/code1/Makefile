# $NEMO/src/nbody/evolve/sellwood/code1
include $(NEMOLIB)/makedefs
#

#  unix ar used, set AR=0 for mac, where this appears not to work
AR=1

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

galaxy:  galaxy$(AR)



galaxy1: galaxy.a $(SRC) $(INC)
	$(FC) $(FFLAGS) -o galaxy galaxy.f galaxy.a

# go figure, on a mac galaxy.a is claimed to be not of the right architecture
# this kludge will create it
galaxy0: $(SRC) $(INC)
	$(FC) $(FFLAGS) -o galaxy galaxy.f lib/*.o


galaxy.a: $(INC)
	(cd lib; make lib)

install:	galaxy
	cp $? $(NEMOBIN)
