#	Makefile for NEMO/src/nbody/evolve/aarseth
include $(NEMOLIB)/makedefs
#

#	Directories to be visited, in this order
DIRS = nbody0 nbody1 nbody2 nbody4 tools # triple

TARNAME = stuff.tar
CHKFILE = Last_update

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
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) bins); done

update:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) update); done

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
