# NEMO/src/Makefile
#	xx-jun-87	Version 1.0	P.J.Teuben
#	 6-jun-88	added lyman and gravsim
# 	 1-jul-89	deleted gravsim, added aarseth
#	   sep-90	added unemo(micro-nemo) in theory
#	   oct-90	new dir.str. for NEMO V2.x

include $(NEMOLIB)/makedefs

DIRS = tools kernel nbody image orbit 

TARDIRS = kernel nbody
TARNAME = kernel.tar

help:
	@echo no help here

dirs:
	-@for i in $(DIRS); do \
	echo $$i; done

cflags:
	-@for i in $(DIRS); do \
	(if [  ! -d $$i ]; then \
           echo -DNO_$$i | tr "[a-z]" "[A-Z]"; \
	fi); done

clean:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) clean); done

cleanlib:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) cleanlib); done

install:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) install); done
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

bins:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) bins); done
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"


update:
	-@for i in $(DIRS); do\
	(cd $$i;$(MAKE) update); done

#	for GNU ld (works for linux)
shared: $(NEMOLIB)/shared shared.$(NEMOHOST)

$(NEMOLIB)/shared:
	mkdir $(NEMOLIB)/shared

clean_shared: shared_clean

shared_clean:
	rm -f $(NEMOLIB)/shared/*

shared.linux:
	(cd $(NEMOLIB); \
		ld -shared -o libnemo.so --whole-archive libnemo.a)

shared.sun5: 
	(cd $(NEMOLIB)/shared; \
		ar x ../libnemo.a; \
		ld -o ../libnemo.so -ztext -G -i *.o)


tar:
	@echo Creating tarfile TARNAME=$(TARNAME) from TARDIRS=$(TARDIRS)
	@echo Starting with kernel/Makefile
	@tar -cf $(TARNAME) kernel/Makefile
	-@for i in $(TARDIRS); do\
	(echo Adding $$i ; tar -rf $(TARNAME) `find $$i \( -name Makefile \
	 -o -name '*.[ch135]' \) -print`);done
	@echo -n Compressing tar file...
	@compress $(TARNAME)
	@echo -n "Done!"

testsuite:
	(cd scripts; ./testsuite -q -i; ./testsuite)
zero:
	@find . -type f -size 0 -print

obj:
	-@for i in $(DIRS);do\
	(find  $$i -type f \( -name \*.o -o -name \*.a -o -name \*.tar \) -print );done

junk:
	-@for i in $(DIRS);do\
	(find  $$i -type f \( \! -name \*.[fhcC] \! -name \*akefile\* \! -name Testfile \! -name .\*.\* \) -print );done

dot:
	-@for i in $(DIRS); do\
	(find $$i -type f -name \.\* -print); done

cc:
	@echo "Your CC compiler is: CC=$(CC)"
	@echo "If it is a shell script, it displays as:"
	@grep -v '#' `which $(CC)`

count:
	wc */*.[fcCh] */*/*.[fcCh] */*/*/*.[Cfch] */*/*/*/*.[Cfch] */*/*/*/*/*.[Cfch] */*/*/*/*/*/*.[Cfch]

count1:
	cat `find . -name \*.c -o -name \*.h -o -name \*.f -o -name \*.C `| wc -l 

