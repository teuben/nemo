MAN8 = *.8
MAN3 = *.3
MAN1 = *.1
MAN5 = *.5
TEX = nemo.tex summary.tex
DVI = nemo.dvi summary.dvi

install: install_man install_tex install_dvi

install_man: install_man8 install_man5 install_man3 install_man1 

install_man8:	$(MAN8)
	cp $? $(NEMO)/man/man8
	@touch install_man8

install_man5:	$(MAN5)
	cp $? $(NEMO)/man/man5
	@touch install_man5

install_man3:	$(MAN3)
	cp $? $(NEMO)/man/man3
	@touch install_man3

install_man1:	$(MAN1)
	cp $? $(NEMO)/man/man1
	@touch install_man1

install_tex:	$(TEX) mkmanual fig*
	cp $? $(NEMO)/text
	@touch install_tex

install_dvi:	$(DVI)
	cp $? $(NEMO)/text
	@touch install_dvi

summary.dvi:  summary.tex
	latex summary.tex

nemo.dvi:  nemo.tex
	@echo "Updating NEMO manual with script 'mkmanual'"
	@echo "This will take a while..."
	mkmanual

clean:
	rm -f *.dvi *.ps *.toc *.aux *.lof *.lot *.idx *.ind *.log
