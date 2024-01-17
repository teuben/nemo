export OPENGL=0# Set this to 1 to enable OpenGL
export SERVER=1# Set this to 1 to enable the visualization web server
include ../rebound/src/Makefile.defs

BINS = reb_viewer reb2s s2reb reb_integrate


# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
CCPROBLEM0=$(CC) -I$(NEMOINC) -I$(NEMOLIB) $(OPT) $(PREDEF) $< -L$(NEMOLIB) -lnemo -lrebound $(LIB)

.PHONY: all librebound


all:   $(BINS)

install: $(BINS)
	mv $(BINS) $(NEMOBIN)

reb_viewer: reb_viewer.c
	$(CCPROBLEM0) -o reb_viewer

reb_integrate: reb_integrate.c
	$(CCPROBLEM0) -o reb_integrate


reb2s: reb2s.c
	$(CCPROBLEM0) -o reb2s

s2reb: s2reb.c
	$(CCPROBLEM0) -o s2reb

clean:
	@echo No clean here yet
