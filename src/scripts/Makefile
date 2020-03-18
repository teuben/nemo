# 	$NEMO/src/scripts/Makefile

include $(NEMOLIB)/makedefs


SCRIPTS = mknemo mkpdoc manlaser nemoman ldso nemo.version \
	fitsedit fitsmerge nds9 shtool nemo.coverage showvar mir2nemo getline pedit \
	wgetc
MAKES = Makefile.lib 

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


scripts: $(SCRIPTS)
	cp $? $(NEMOBIN)
	cp Makefile.lib $(NEMOLIB)/Makefile
