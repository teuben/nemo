# 	$NEMO/src/scripts/Makefile

include $(NEMOLIB)/makedefs


SCRIPTS = mknemo mkpdoc manlaser nemoman nemobench ldso nemo.version \
	fitsedit fitsmerge nds9 shtool nemo.coverage showvar mir2nemo getline pedit \
	wgetc yapp_query nemoplot nemopars nemo \
        sbatch_nemo.sh nemo_functions.sh progress_bar.sh \
	python/tabplot.py python/fitsplot.py
MAKES = Makefile.lib 

help:
	@echo special NEMO stuff
	@echo No more specific help here

####		install by NEMO (done by top level directly)

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


checkpars:
	(cd ../../; src/scripts/checkpars.py -v)  > checkpars.log
	@grep ^'Files read'         checkpars.log
	@grep ^'Bad files found'    checkpars.log
