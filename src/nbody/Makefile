#	Makefile for NEMO/src/nbody
include $(NEMOLIB)/makedefs

#	Directories to be visited, in this order
DIRS = cores cores/bodysub io io_nemo tools init evolve trans reduc \
	fits kepler image xyz

TARNAME = stuff.tar
CHKFILE = Last_update

#.IGNORE:  glnemo

help:
	@echo NEMO V2.2 NEMO/src/nbody

clean:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) clean); done

cleanlib:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) cleanlib); done

install:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) install); done

bins:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) bins); done
	-(cd glnemo; qmake; make install)

update:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) update); done

tar:
	@echo Creating tarfile TARNAME=$(TARNAME) from DIRS=$(DIRS)
	@echo -n "Hang on..."
	@tar cf $(TARNAME) Makefile
	-@for i in $(DIRS); do\
	(tar rf $(TARNAME) `find $$i \( -name Makefile \
	 -o -name '*.[ch135]' \) -print`);done
	@echo  "and compressing $(TARNAME) ...."
	@compress $(TARNAME)
	@echo "  Done!"

newtar:
	@tar cvf $(TARNAME) `find . -newer $(CHKFILE) -type f -print`	

newdir:
	@find . -newer $(CHKFILE) -type f -exec ls -l '{}' \;

glnemo:
	@(cd glnemo; qmake; make)
