#	Makefile for NEMO/src/nbody/evolve/sellwood
#       this is an old version (V1.3) of the galaxy code.
#       See comments in $NEMO/usr/sellwood how to install a newer release.
include $(NEMOLIB)/makedefs
#

#	Directories to be visited, in this order
DIRS = code1 tools


clean:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) clean); done

cleanlib:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) cleanlib); done

install:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) install); done

bins:
	@echo No bins in sellwood

update:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) update); done

