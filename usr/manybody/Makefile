

CVS1 = partiview  
CVS2 = staratlas  starcluster  stardata  starlab

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
