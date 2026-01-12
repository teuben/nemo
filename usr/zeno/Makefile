# 
URL0 = http://legacy.ifa.hawaii.edu/faculty/barnes/
URL1 = https://github.com/joshuabarnes/zeno/
URL2 = https://github.com/teuben/zeno

TIME = /usr/bin/time

UBUNTU = libgsl-dev openmpi-bin freeglut3-dev libopenmpi-dev


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
	git clone $(URL1) zeno_jeb

#  my current repo has a branch for linux
zeno_pjt:
	git clone $(URL2) zeno_pjt
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

#  install for ubuntu
ubuntu:
	sudo apt install $(UBUNTU)

# IFA has changed their access
WFLAGS = --no-check-certificate


idkit:  identikit2

identikit2:
	wget $(WFLAGS) $(URL0)/research/identikit/identikit2.tar
	mkdir identikit2
	(cd identikit2; tar xf ../identikit2.tar)

testdata.txt:
	wget $(WFLAGS) $(URL0)/treecode_old/testdata.txt

code1:
	(cd code_c; make clean all)

bench.ascii:
	mkplummer - 128 | stoa - - > bench.ascii

bench1: bench.ascii
	code_c/code < bench.in > bench.log   2 > bench.err

# benchmark noted in bhtree(1NEMO)
bench2:
	@rm -f p1k p1k.out
	plumsnap p1k 1024
	$(TIME) treecode    p1k p1k.out dtime=1/64 theta=0.75 usequad=false tstop=10 > /dev/null
	$(TIME) treecode    p1k p1k.out dtime=1/64 theta=0.75 usequad=true  tstop=10 > /dev/null
	$(TIME) treecode_dp p1k p1k.out dtime=1/64 theta=0.75 usequad=false tstop=10 > /dev/null
	$(TIME) treecode_dp p1k p1k.out dtime=1/64 theta=0.75 usequad=true  tstop=10 > /dev/null
