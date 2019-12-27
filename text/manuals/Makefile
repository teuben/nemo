# 	Makefile for Miriad LaTex manuals - only for the experts
#	xx-feb-91   pjt      Created
#       21-jul-93   pjt      Added 'src' verbatim includes concept
#        6-aug-95   pjt      Added html

F   = nemo
TEX = nemo.tex latest.tex idx.tex
STY = nemo.sty troffman.sty html.sty
INC = ack.inc aiface.inc bench.inc coordsys.inc dirs.inc disclaimer.inc examples.inc \
      filestr.inc glossary.inc graphics.inc gsl.inc history.inc iface.inc \
      image.inc \
      install.inc intro.inc miriad.inc mirtool.inc \
      orbit.inc potential.inc potname.inc potctex.inc \
      progr.inc readme.inc refs.inc setup.inc shell.inc \
      trouble.inc convent.inc \
      system.inc bodies.inc \
      psinsert.inc
SRC = hello.src makefile.src snapfirst.src testf2c.src nemo.ind 
PS  = encounter1.ps plummer.ps
TMP = 
EXT =
SUB = $(INC) $(STY) $(TEX) $(EXT) $(SRC) $(PS) Makefile README


# for submission into the NEMO tree
SUBMIT = nemo.submit
SFLAGS =
DIR    = $(NEMO)/text/manuals

help:
	@echo "         --- be sure to read the README file in this dir"
	@echo "Make $(F) manual etc. using latex"
	@echo "================================="
	@echo "Macro names:  "
	@echo "  F, TEX, STY, INC, EXT, SUB, DIR, SUBMIT"
	@echo "Targets: "
	@echo "  help		this help"
	@echo "  cleanall       clean all .dvi, .aux etc. tex related thingos"
	@echo "  clean [F=$(F)] like cleanall, but only for $(F).* files"
	@echo "  $(F)           make dvi file for manual F=$(F)"
	@echo "  submit         submit text to DIR=$(DIR) using SUBMIT=$(SUBMIT)"
	@echo "  "
	@echo "The used INC= files in this dir are:"
	@echo "  $(INC)"
	@echo "The used TEX= files in this dir are:"
	@echo "  $(TEX)"
	@echo "The used STY= files in this dir are:"
	@echo "  $(STY)"
	@echo "The used EXT= files in this dir are:"
	@echo "  $(EXT)"
	@echo "The used TMP= files in this dir are:"
	@echo "  $(TMP)"



cleanall:
	@echo "Cleaning all (La)TeX related documents of all files.....";sleep 2
	rm -f *.log *.aux *.toc *.dvi *.lot *.lof *.ps *.ind *.idx *.tmp
	rm -rf tmp

clean:
	@echo "Cleaning all (La)TeX related documents for F=$(F)....";sleep 2
	rm -f $(F).log $(F).aux $(F).toc $(F).dvi $(F).lot $(F).lof \
            $(F).ps $(F).ind $(F).idx
	@echo "The .tmp files have not been deleted. Only in 'make cleanall'."

$(F): $(F).dvi $(F).ps

#   to do this right, perhaps check rubber, a latex wrapper
#   http://www.pps.jussieu.fr/~beffara/soft/rubber/
$(F).dvi: $(TEX) $(INC) $(STY)
	@echo "Running latex 3x - will take a while [did you run 'make .prog'?]"
	-latex $(F) 
	-latex $(F)
	-makeindex < $(F).idx > $(F).ind
	-latex $(F)
	-latex $(F)

$(F).ps:  $(F).dvi
	dvips $(F) -o
	@echo Also making a 2-page-on-1 version
	psnup -2 nemo.ps > nemo2.ps


tar:
	tar cvf $(F).tar $(SUB)

submit:  .submit

.submit:    $(SUB)
	$(SUBMIT) $(SFLAGS) $? $(DIR)
	@touch .submit

diff:
	-@for i in $(SUB); do \
	echo ">> $$i <<";diff $$i $(DIR); done
    
cmp:
	-@for i in $(SUB); do \
	cmp $$i $(DIR)/$$i; done

TEST = test.inc
test:
	@echo "Usage: make test TEST=$(TEST)"
	@echo "where TEST is a latex include file,  normally w/ extension .inc"
	echo $(TEST) | latex  latest
	makeindex <latest.idx > latest.ind
	echo $(TEST) | latex  latest

html:
	@echo Creating html manual using latex2html
	@echo You should have done a full latex run before
	-latex2html \
		-address "(c) Peter Teuben" \
		-up_url http://www.astro.umd.edu/nemo/ \
		-up_title "NEMO Home Page" \
		-show_section_numbers \
		-html_version 3.0 \
		nemo.tex


# 		-ps_images \

pdf:
	pdflatex nemo
	pdflatex nemo
	pdflatex nemo
