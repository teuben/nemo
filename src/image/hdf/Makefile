#	Makefile for NEMO/src/image/hdf
include $(NEMOLIB)/makedefs

#	Directories to be visited, in this order
DIRS = cmhog


help:
	@echo Nemo makefile in `pwd`

clean:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) clean); done

cleanlib:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) cleanlib); done

install:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) install); done

update:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) update); done

