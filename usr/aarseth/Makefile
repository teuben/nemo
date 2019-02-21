#  Aarseth, and related code, not formally downloaded

include $(NEMOLIB)/makedefs

URL1 = ftp://ftp.ast.cam.ac.uk/pub/sverre/nbody6/nbody6.tar.gz
URL2 = https://github.com/nbodyx/Nbody6ppGPU

help:
	@echo "install:         install NEMOBIN/gravidy"
	@echo "nbody6:          grab the code from URL1=$(URL1)"
	@echo "nbody6++:	grab the code from URL2=$(URL2)"

nbody6:
	wget $(URL1)
	(mkdir nbody6 ; cd nbody6; tar zxf ../nbody6.tar.gz)

# just an example on a simple laptop
install_nbody6:
	(cd nbody6/GPU2 ;\
	make sse ;\
	cp run/nbody7.sse $(NEMOBIN)/nbody6)



nbody6++: Nbody6ppGPU

Nbody6ppGPU:
	git clone $(URL2)

clean:
	rm -rf gravidy

# just an example on a simple laptop
install_nbody6++:
	(cd Nbody6ppGPU ;\
	./configure --disable-gpu --disable-mpi ;\
	make ;\
	cp build/nbody6++.avx $(NEMOBIN)/nbody6++)

install_nbody6:

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
