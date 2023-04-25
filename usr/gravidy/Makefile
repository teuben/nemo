include $(NEMOLIB)/makedefs

URL1 =  http://github.com/cmaureir/gravidy.git
URL2 =  http://github.com/chongchonghe/gravidy.git
URL3 =  http://github.com/teuben/gravidy.git
URL4 =  https://gitlab.com/cmaureir/gravidy-view

.PHONY : help
## help:         This Help
help : Makefile
	@sed -n 's/^##//p' $<

## url:          show the convenience URL's we have defined
url:
	@echo "URL1:   $(URL1)"
	@echo "URL2:   $(URL2)"
	@echo "URL3:   $(URL3)"
	@echo "URL4:   $(URL4)"

## gravidy:      get the github repo from URL1
gravidy:
	git clone $(URL1)

## clean:        distclean the gravidy source tree
clean:
	(cd gravidy/src/; make distclean)

## cleanall:     remove data creates by local tests
cleanall:
	rm -rf p128 run3 run4 run5

## install:      install the gravidy-cpu into $NEMOBIN
install: gravidy
	(cd gravidy/src ; make clean pn $(NEMO_BFLAGS); cp gravidy-cpu $(NEMOBIN)/gravidy)

## gravidy-view: get the 3D gravidy git repo. Install at your leisure
gravidy-view:
	git clone $(URL4)

## ubuntu:       example ubuntu-style install target to get boost
ubuntu:
	sudo apt install libboost-all-dev


NEMO_BFLAGS = "BOOSTFLAGS =-I$(NEMO)/opt/boost_1_69_0 -L$(NEMO)/opt/boost_1_69_0/stage/lib -lboost_program_options"
NEMO_NVCC   = /test/usr/local/cuda/bin/nvcc
## install_gpu:  example gravidy-gpu
install_gpu: gravidy
	(cd gravidy/src ; \
	make clean gpu NVCC=$(NEMO_NVCC) BOOSTFLAGS = "$(NEMO_BFLAGS)"; \
	cp gravidy-cpu $(NEMOBIN)/gravidy)	



# example runs

NBODY = 128

## p128:         benchmark plummer sphere osf 128 particles for regression
p128:
	mkplummer p128 $(NBODY) seed=128

## run3:         gravidy run
run3:   p128
	mkdir run3
	snapprint p128 i,m,x,y,z,vx,vy,vz format=%.15f > run3/in
	gravidy --input run3/in --output run3/run3 --interval 0.1 --time 10 --eta 0.01 --softening 0.05
	cat run3/run3.out.snapshot_0* | sed 's/# Time://' |\
                tabtos - - time skip,m,pos,vel nbody=128 |\
                csf - run3/OUT3.snap SnapShot
## run3:         rungravidy run
run4:	p128
	rungravidy p128 run4 eta=0.01 deltat=0.1 tcrit=10.0 eps=0.05

## run3:         hackcode1 run
run5:	p128
	mkdir run5
	hackcode1 p128 run5/OUT3.snap tstop=10 eps=0.05 freqout=10 freq=100

## cmp34:        compare run3 and run4
cmp34:  run3 run4
	snapcmp run3/OUT3.snap run4/OUT3.snap

## cmp45:        compare run4 and run5
cmp45:  run4 run5
	snapcmp run4/OUT3.snap run5/OUT3.snap
