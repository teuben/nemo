
EXPORT = tkrun.c tkrun.1 tkrun.lsm Makefile \
	testscript testscript0.tk redo \
	sample1 \
	CHANGES COPYING README TODO \
	adass97.ps.gz
VERSION = 0.2e
FTPDIR = apus.astro.umd.edu:/home/ftp/pub/nemo

CFLAGS = -g

help:
	@echo There is no makefile help here

clean:
	rm -f core *.o tkrun testscript.tk  scratchrun.tk


tkrun: tkrun.c
	$(CC) $(CFLAGS) -o tkrun tkrun.c


export: $(EXPORT)
	tar cf tkrun-$(VERSION).tar $(EXPORT)
	-mkdir tkrun-$(VERSION)
	(cd tkrun-$(VERSION) ; tar xf ../tkrun-$(VERSION).tar)
	tar cf tkrun-$(VERSION).tar tkrun-$(VERSION)
	gzip tkrun-$(VERSION).tar 
	cp tkrun-$(VERSION).tar.gz $(FTPDIR)

export_diff:
	-@for i in $(EXPORT); do\
	(echo "###: $$i ";diff $$i tkrun-$(VERSION)); done


diff: tkrun
	@echo Generating a new testscript.tk
	@tkrun -x testscript > /dev/null 2>&1
	@echo And running diff:
	diff testscript0.tk testscript.tk 

update:
	cp testscript.tk  testscript0.tk 


dist:
	cvs -q export -D tomorrow -d tkrun-$(VERSION) tkrun
	tar zcf tkrun-$(VERSION).tar.gz tkrun-$(VERSION)
	rm -rf tkrun-$(VERSION)
	scp tkrun-$(VERSION).tar.gz $(FTPDIR)

