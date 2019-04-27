#  Aarseth, and related code, with download instructions to install

include $(NEMOLIB)/makedefs

URL1 = ftp://ftp.ast.cam.ac.uk/pub/sverre/nbody6/nbody6.tar.gz
URL2 = https://github.com/nbodyx/Nbody6ppGPU

help:
	@echo "install:           install NEMOBIN/gravidy"
	@echo "nbody6:            grab the code from URL1=$(URL1)"
	@echo "nbody6++:          grab the code from URL2=$(URL2)"
	@echo "install_nbody6:    install in $(NEMOBIN)"
	@echo "install_nbody6++:  install in $(NEMOBIN)"

nbody6:
	wget $(URL1)
	(mkdir nbody6 ; cd nbody6; tar zxf ../nbody6.tar.gz)

# just an example on a simple laptop
# note the order of building appears to matter !!!
install_nbody6:  nbody6
	(cd nbody6/GPU2 ;  make sse ; cp run/nbody7.sse $(NEMOBIN))
	(cd nbody6/Ncode ; make ;     cp nbody6 $(NEMOBIN))


nbody6++: Nbody6ppGPU

Nbody6ppGPU:
	git clone $(URL2)


# just an example on a simple laptop
install_nbody6++:  Nbody6ppGPU
	(cd Nbody6ppGPU ;\
	./configure --disable-gpu --disable-mpi ;\
	make ;\
	cp build/nbody6++.avx $(NEMOBIN)/nbody6++)

# example runs

NBODY = 1024
TCRIT = 1.0

p.snap:
	mkplummer p.snap $(NBODY) seed=123

# nbody6++
run1:   p.snap
	time runbody6 p.snap run1 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05

# nbody6 single processor
run2:   p.snap
	time runbody6 p.snap run2 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=1 exe=nbody6 

# nbody6 SSE multiprocessor
run3:   p.snap
	time runbody6 p.snap run3 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=1 exe=nbody7.sse

run4:	p.snap
	time runbody4 p.snap run4 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 

#       for kicks, a softened gravity run
run0:
	time rungravidy p.snap run0 tcrit=$(TCRIT) deltat=0.1 eps=0.01

