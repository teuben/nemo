
CVS1 = partiview  
CVS2 = staratlas  starcluster  stardata  starlab

CVSROOT = :pserver:anonymous@cvs.astro.umd.edu:/home/cvsroot
COPT = -Q -d $(CVSROOT)

echo:
	@echo This directory is a placeholder for the manybody tools,
	@echo sans NEMO of course. 

cvs:	cvs1 cvs2

cvs1:
	-@for i in $(CVS1); do \
	(echo $$i; cvs -Q co $$i); done

cvs2:
	-@for i in $(CVS2); do \
	(echo $$i; cvs -Q co -d $$i manybody/$$i); done

clean:
	rm -rf $(CVS1) $(CVS2)

junk:
	cvs -Q -d $(ACVSROOT) co -d nemo_cvs nemo
	date > nemo_cvs/VERSION_cvs
	echo "$(ACVSROOT) A" > nemo_cvs/.cvspass
	cp nemo_cvs/VERSION_cvs /tmp
	gtar zcf nemo_cvs.tar.gz nemo_cvs
	rm -rf nemo_cvs
	scp nemo_cvs.tar.gz /tmp/VERSION_cvs apus:/home/ftp/pub/nemo

cvstar: cvstar1 cvstar2 cvstar3 cvstar4

cvstar1: 
	rm -rf nemo
	cvs $(COPT) co  nemo;  (cd nemo; mkdir local; cd local; cvs $(COPT) co pgplot)
	gtar zcf nemo.tar.gz nemo

cvstar2: 
	rm -rf partiview
	cvs $(COPT) co  partiview; gtar zcf partiview.tar.gz partiview

cvstar3: 
	rm -rf starlab
	cvs -Q co  starlab;   gtar zcf starlab.tar.gz starlab

cvstar4: 
	rm -rf starcluster
	cvs -Q co  starcluster;   gtar zcf starcluster.tar.gz starcluster



	
