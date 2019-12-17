# $NEMO/src/nbody/evolve/sellwood/code1
# now produces the old galaxy13
include $(NEMOLIB)/makedefs
#

#  unix ar used, set AR=0 for mac, where this appears not to work
#  so check for Darwin
SYS = `uname -s`

INC = admin.h  rjinclude.h
SRC = galaxy.f

help:
	@echo GALAXY program, V1.3 written by Jerry Sellwood, for SYS=$(SYS)

#   on 64bit machines, -O will coredump, -O2 is ok, adding -g also
#   we'll inherit FFLAGS from the makedefs file instead
#FFLAGS = -O2

all:	galaxy.a galaxy13

clean:
	rm -f $(LIB) galaxy *.o *.a core galaxy.res galaxy.fin galaxy.lis
	(cd lib; make clean)




# go figure, on a mac galaxy.a is claimed to be not of the right architecture
# this kludge will create it. Will 'libtool -static' work instead of 'ar' on mac ?


galaxy13: lib0 $(SRC) $(INC) 
	$(FC) $(FFLAGS) -o galaxy13 galaxy.f lib/*.o

lib0:
	(cd lib; make lib)


# although this seems standard unix, doesn't work on a mac
galaxy1: galaxy.a $(SRC) $(INC)
	$(FC) $(FFLAGS) -o galaxy13 galaxy.f galaxy.a


galaxy.a: $(INC)
	(cd lib; make lib)

install:	galaxy13
	cp $? $(NEMOBIN)
