include $(NEMOLIB)/makedefs

help:
	@echo "install:         install NEMOBIN/gravidy"
	@echo "gravidy:         grab the code"
	@echo "ubuntu:          install BOOST (ubuntu)"
	@echo "gravidy-view:    install the 3D viewer"

gravidy:
	git clone http://github.com/cmaureir/gravidy.git

clean:
	rm -rf gravidy

install: gravidy
	(cd gravidy/src ; make clean pn $(NEMO_BFLAGS); cp gravidy-cpu $(NEMOBIN)/gravidy)

gravidy-view:
	git clone https://gitlab.com/cmaureir/gravidy-view

ubuntu:
	sudo apt install libboost-all-dev


NEMO_BFLAGS = "BOOSTFLAGS =-I$(NEMO)/opt/boost_1_69_0 -L$(NEMO)/opt/boost_1_69_0/stage/lib -lboost_program_options"
NEMO_NVCC   = /test/usr/local/cuda/bin/nvcc
install_gpu: gravidy
	(cd gravidy/src ; \
	make clean gpu NVCC=$(NEMO_NVCC) BOOSTFLAGS = "$(NEMO_BFLAGS)"; \
	cp gravidy-cpu $(NEMOBIN)/gravidy)	



# example runs

NBODY = 128

p128:
	mkplummer p128 $(NBODY) seed=128


run3:   p128
	mkdir run3
	snapprint p128 i,m,x,y,z,vx,vy,vz format=%.15f > run3/in
	gravidy --input run3/in --output run3/run3 --interval 0.1 --time 10 --eta 0.01 --softening 0.05
	cat run3/run3.out.snapshot_0* | sed 's/# Time://' |\
                tabtos - - time skip,m,pos,vel nbody=128 |\
                csf - run3/OUT3.snap SnapShot

run4:	p128
	rungravidy p128 run4 eta=0.01 deltat=0.1 tcrit=10.0 eps=0.05

cmp34:  run3 run4
	snapcmp run3/OUT3.snap run4/OUT3.snap
