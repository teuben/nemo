
#

zeno:  zeno_jeb zeno_pjt
	@echo 'You will need to create a symlink zeno -> zeno_XYZ'
	@echo Pick from:
	@echo ln -s zeno_jeb zeno
	@echo ln -s zeno_pjt zeno

jeb:
	rm -f zeno; ln -s zeno_jeb zeno

pjt:
	rm -f zeno; ln -s zeno_pjt zeno

#  original one
zeno_jeb:
	git clone https://github.com/joshuabarnes/zeno/ zeno_jeb

#  my current repo has a branch for linux
zeno_pjt:
	git clone https://github.com/teuben/zeno zeno_pjt
	(cd zeno_pjt; git co teuben1)

#  install
install:
	(source zeno_start.sh ; cd zeno; make -f Zeno install_all)
	@echo 'Consider "make zbin" to isolate namesake ZENO bins from NEMO bins:'
	@echo ZBIN=$(ZBIN)

#  move ZENO bin's to bin/zeno if their name clashes with NEMO
ZBIN = tsf csf rsf \
	snapkey snaprect snaprotate snapsort snapsplit snapstack snapstat
zbin:
	mkdir -p zeno/bin/zeno
	(cd zeno/bin; mv $(ZBIN) zeno)

#  optional install of gsl in $NEMO/opt
gsl:
	mknemo gsl

# IFA has changed their access
WFLAGS = --no-check-certificate
URL1 = http://legacy.ifa.hawaii.edu/faculty/barnes/

idkit:  identikit2

identikit2:
	wget $(WFLAGS) $(URL1)/research/identikit/identikit2.tar
	mkdir identikit2
	(cd identikit2; tar xf ../identikit2.tar)

testdata.txt:
	wget $(WFLAGS) $(URL1)/treecode_old/testdata.txt

code1:
	(cd code_c; make clean all)

bench.ascii:
	mkplummer - 128 | stoa - - > bench.ascii

bench1: bench.ascii
	code_c/code < bench.in > bench.log   2 > bench.err
