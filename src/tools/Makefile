#	Makefile for NEMO/src/kernel
include $(NEMOLIB)/makedefs

#	Directories to be visited, in this order
#DIRS = movie movietool
DIRS = misc tkrun


TARNAME = stuff.tar
CHKFILE = Last_update

help:
	@echo Nemo makefile in `pwd`
	@echo Targets:
	@echo "   clean, install, update"
	@echo "   tar     source in TARNAME=$(TARNAME) from DIRS=$(DIRS)"
	@echo "   dir     list source for tar"
	@echo "   newtar  new stuff in TARNAME=$(TARNAME) from
	@echo "   newdir  list new stuff"

clean:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) clean); done

cleanlib:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) cleanlib); done

install:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) install); done

bins:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) bins); done

update:
	-@for i in $(DIRS); do\
	(cd $$i;echo "### Dir: `pwd`";$(MAKE) $(MFLAGS) update); done

tar:
	@echo Creating tarfile TARNAME=$(TARNAME) from DIRS=$(DIRS)
	@echo -n "Hang on..."
	@tar cf $(TARNAME) Makefile ../NEMO-boot
	-@for i in $(DIRS); do\
	(tar rf $(TARNAME) `find $$i \( -name Makefile \
	 -o -name '*.[ch135]' \) -print`);done
	@echo  "and compressing $(TARNAME) ...."
	@compress $(TARNAME)
	@echo "  Done!"
dir:
	@ls Makefile
	-@for i in $(DIRS); do\
	(find $$i \( -name Makefile -o -name '*.[ch135]' \) -print);done

newtar:
	@tar cvf $(TARNAME) `find . -newer $(CHKFILE) -type f -print`	

newdir:
	@find . -newer $(CHKFILE) -type f -exec ls -l '{}' \;

obj:
	find . -type f -name \*.o -print



# special targets, to make life easy

tr2latex: .tr2latex

.tr2latex:
	@echo Making tr2latex 
	(cd tr2latex; make -f Makefile.unix; mv tr2latex $(NEMOBIN); make -f Makefile.unix clean)
