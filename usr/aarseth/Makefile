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


# just an example on a simple laptop
install_nbody6++:
	(cd Nbody6ppGPU ;\
	./configure --disable-gpu --disable-mpi ;\
	make ;\
	cp build/nbody6++.avx $(NEMOBIN)/nbody6++)

install_nbody6:

# example runs

NBODY = 1024

p.snap:
	mkplummer p.snap $(NBODY) seed=123


run1:   p.snap
	runbody6 p.snap run1 tcrit=1 deltat=0.1

run2:   p.snap
	runbody6 p.snap run2 tcrit=1 deltat=0.1 exe=nbody6 nbody6=0
