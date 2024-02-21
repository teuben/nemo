# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# PJT keeps a private version because he doesn't trust the cloud (and has been proven right in July 2021)
pjt: lastmod html
	rsync -av _build/html  ssh.astro.umd.edu:public_html/nemo/readthedocs
	@echo https://www.astro.umd.edu/~teuben/nemo/readthedocs/html

# PJT2 in case latexmk was not installed, the "make latexpdf" may fail, and manually may work here
pjt2: lastmod html
	make latexpdf; cd _build/latex/; pdflatex nemo.tex
	@echo xdg-open _build/latex/nemo.pdf

lastmod:
	echo "Last modified `date`" > $(SOURCEDIR)/lastmod.rst
#	echo "`git rev-list --count HEAD` and Revision `git rev-parse master`" >>  $(SOURCEDIR)/lastmod.rst
