# 	$NEMO/src/scripts/Makefile

include $(NEMOLIB)/makedefs

#	default extension for cc and make scripts (not used anymore)
OS=
#
#CFLAGS = -O
#LDFLAGS = -O

SCRIPTS = mknemo mkpdoc manlaser nemoman ldso nemo.version \
	fitsedit nds9 shtool nemo.coverage showvar mir2nemo getline
MAKES = Makefile.lib 
FAKES = cc make f77 ranlib

help:
	@echo special NEMO stuff
	@echo No more specific help here

####		install by NEMO (done by top level directly)
#
# CC is somewhat of an exception: it should be installed before
# everything else by a seperate call from the top-level Makefile
#
#nemo_startup: cc make
#	cp $? $(NEMOBIN)
#
####
# mycc.c, which should be still in this directory, has become obsolete
# when using the shell script, Ultrix does not support -L flag, causing
# Makefiles to explicitely name the libraries.....

nemo_lib:
	cp Makefile.lib $(NEMOLIB)/Makefile

nemo_bin: $(SCRIPTS)
	cp $? $(NEMOBIN)
	cp BADPROGRAM $(NEMOBIN)
	cp BADPROGRAM $(NEMO)/bin
	rm -f *.o *.a core make cc


nemo_src:

clean:
	rm -f *.o *.a core 
####				install by USER (nemo)

install: install_lib install_bin
	@echo utilities installed

install_lib: $(MAKES)
	cp $? $(NEMOLIB)

install_bin: $(SCRIPTS)
	cp $? $(NEMOBIN)

tidy:
	rm -f *.o *.BAK *.CKP core

mycc:	mycc.c
	$(CC) $(CFLAGS) -o cc mycc.c


fakes:	$(FAKES)

#
#  special cc and make scripts; called by root Makefile during install
#

cc:
	@-if [ -f cc.d/cc.$(OS) ]; \
	then \
	  echo "Copying cc.d/cc.$(OS) script as cc replacement to NEMOBIN"; \
	  if [ ! -f $(NEMOBIN)/cc ]; \
	  then \
	    cp cc.d/cc.$(OS) $(NEMOBIN)/cc; \
	  else \
	    echo Will not overwrite existing $(NEMOBIN)/cc; \
	  fi; \
	else \
	  echo "No cc script available for OS=$(OS) - skipped"; \
	  echo "Standard 'cc' will hence be used:"; \
	  which cc; \
	fi

f77:
	@-if test -f f77.d/f77.$(OS); \
	then \
	  echo "Copying f77.d/f77.$(OS) script as f77 replacement to NEMOBIN"; \
	  if [ ! -f $(NEMOBIN)/f77 ]; \
	  then \
	    cp f77.d/f77.$(OS) $(NEMOBIN)/f77; \
	  else \
	    echo Will not overwrite existing $(NEMOBIN)/f77; \
	  fi; \
	else \
	  echo "No f77 script available for OS=$(OS) - skipped"; \
	  echo "Standard f77 compiler will hence be used: "; \
	  which f77; \
	fi

ranlib:
	@-if test -f ranlib.d/ranlib.$(OS); \
	then \
	  echo "Copying ranlib.d/ranlib.$(OS) script as ranlib replacement to NEMOBIN"; \
	  if [ ! -f $(NEMOBIN)/ranlib ]; \
	  then \
	    cp ranlib.d/ranlib.$(OS) $(NEMOBIN)/ranlib; \
	  else \
	    echo Will not overwrite existing $(NEMOBIN)/ranlib; \
	  fi; \
	else \
	  echo "No ranlib script available for OS=$(OS) - skipped"; \
	  echo "Standard ranlib will hence be used: "; \
	  which ranlib; \
	fi

make:
	@-if test -f make.d/make.$(OS); \
	then \
	  echo "Copying make.$(OS) script as make replacement to NEMOBIN"; \
	  if [ ! -f $(NEMOBIN)/make ]; \
	  then \
	    cp make.d/make.$(OS) $(NEMOBIN)/make; \
	  else \
	    echo Will not overwrite existing $(NEMOBIN)/make; \
	  fi; \
	else \
	  echo "No make script available for OS=$(OS) - skipped"; \
	  echo "Standard /bin/make will hence be used: "; \
	  which make; \
	fi

scripts: $(SCRIPTS)
	cp $? $(NEMOBIN)
	cp Makefile.lib $(NEMOLIB)/Makefile
