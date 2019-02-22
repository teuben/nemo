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
# note the order of building appears to matter !!!
install_nbody6:
	(cd nbody6/GPU2 ;  make sse ; cp run/nbody7.sse $(NEMOBIN)/nbody6.sse)
	(cd nbody6/Ncode ; make ;     cp nbody6 $(NEMOBIN)/nbody6)


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
TCRIT = 1.0

p.snap:
	mkplummer p.snap $(NBODY) seed=123


run1:   p.snap
	time runbody6 p.snap run1 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05

run2:   p.snap
	time runbody6 p.snap run2 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=0 exe=nbody6 

run3:   p.snap
	time runbody6 p.snap run3 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=0 exe=nbody6.sse
