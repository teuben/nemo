include $(NEMOLIB)/makedefs

URL1 =  http://github.com/cmaureir/gravidy.git
URL2 =  http://github.com/chongchonghe/gravidy.git
URL3 =  -b teuben1 http://github.com/teuben/gravidy.git
URL4 =  https://gitlab.com/cmaureir/gravidy-view
URL5 =  http://gitlab.com/cmaureir/gravidy-input.git

URL  = $(URL3)

TIME = /usr/bin/time

.PHONY : help
## help:          This Help
help : Makefile
	@sed -n 's/^##//p' $<

## url:           show the convenience URL's we have defined
url:
	@echo "URL1:   $(URL1)"
	@echo "URL2:   $(URL2)"
	@echo "URL3:   $(URL3)"
	@echo "URL4:   $(URL4)"
	@echo "URL5:   $(URL5)"
	@echo "Be aware some may need a specific branch"

## gravidy:       get the github repo from URL2 for now where we do 2023 development
gravidy:
	git clone $(URL)

## gravidy-input: get sample input files
gravidy-input:
	git clone $(URL5)

## gravidy-view:  get the 3D gravidy git repo. Install at your leisure
gravidy-view:
	git clone $(URL4)

## clean:         distclean the gravidy source tree
clean:
	(cd gravidy/src/; make distclean)

## cleanall:      remove data creates by local tests
cleanall:
	rm -rf p128 run3 run4 run5

## ubuntu:        example ubuntu-style install target to get boost and openmpi
ubuntu:
	sudo apt install libboost-all-dev openmpi-bin


## install:       install rungravidy and the default gravidy-cpu
install: install_cpu
	mknemo rungravidy

## install_cpu:   example gravidy-cpu (uses CXX=)
install_cpu: gravidy
	(cd gravidy/src ;\
	make clean cpu CXX=$(CXX); \
	cp gravidy-cpu $(NEMOBIN)/gravidy; cp gravidy-cpu $(NEMOBIN))


NEMO_BFLAGS = "BOOSTFLAGS =-I$(NEMO)/opt/boost_1_69_0 -L$(NEMO)/opt/boost_1_69_0/stage/lib -lboost_program_options"
NEMO_NVCC   = /test/usr/local/cuda/bin/nvcc
## install_gpu:   example gravidy-gpu (uses NEMO_BFLAGS and NEMO_NVCC)
install_gpu: gravidy
	(cd gravidy/src ; \
	make clean gpu NVCC=$(NEMO_NVCC) BOOSTFLAGS = "$(NEMO_BFLAGS)"; \
	cp gravidy-cpu $(NEMOBIN)/gravidy; cp gravidy-cpu $(NEMOBIN))

## install_mpi:   example gravidy-mpi
install_mpi: gravidy
	(cd gravidy/src ; \
	make clean mpi CXX=$(CXX) ; \
	cp gravidy-mpi $(NEMOBIN)/gravidy; cp gravidy-mpi $(NEMOBIN))



# example runs

NBODY = 128

## p128:          benchmark plummer sphere of 128 particles for regression
p128:
	mkplummer p128 $(NBODY) seed=128

EXE = gravidy

## run3:          gravidy run
run3:   p128
	mkdir run3
	snapprint p128 i,m,x,y,z,vx,vy,vz format=%.15f > run3/in
	$(EXE) --input run3/in --output run3/run3 --interval 0.1 --time 10 --eta 0.01 --softening 0.05
	cat run3/run3.out.snapshot_0* | sed 's/# Time://' |\
                tabtos - - time skip,m,pos,vel nbody=128 |\
                csf - run3/OUT3.snap SnapShot
	bsf run3/OUT3.snap "0.00679654 0.736003 -8.00616 10 91494"

## run4:          rungravidy run
run4:	p128
	rungravidy p128 run4 eta=0.01 deltat=0.1 tcrit=10.0 eps=0.05

## run5:          hackcode1 run
run5:	p128
	mkdir run5
	hackcode1 p128 run5/OUT3.snap tstop=10 eps=0.05 freqout=10 freq=100

## cmp34:         compare run3 and run4
cmp34:  run3 run4
	snapcmp run3/OUT3.snap run4/OUT3.snap

## cmp45:         compare run4 and run5
cmp45:  run4 run5
	snapcmp run4/OUT3.snap run5/OUT3.snap

## example1:      quick example from the manual
example1: gravidy-input
	$(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 1 

## bench1:        single core performance
bench1: gravidy-input
	OMP_NUM_THREADS=1 $(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 10 > /dev/null

## bench2:        multi core performance
bench2: gravidy-input
	OMP_NUM_THREADS=1 $(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 10 > /dev/null
	OMP_NUM_THREADS=2 $(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 10 > /dev/null
	OMP_NUM_THREADS=4 $(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 10 > /dev/null
	OMP_NUM_THREADS=8 $(TIME) $(EXE) -i gravidy-input/04-nbody-p1024_m1.in -p -t 10 > /dev/null
