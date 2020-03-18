#  Aarseth, and related code, with download instructions to install

include $(NEMOLIB)/makedefs

WGET = wget

URL1 = ftp://ftp.ast.cam.ac.uk/pub/sverre/nbody6/nbody6.tar.gz
URL2 = https://github.com/nbodyx/Nbody6ppGPU

help:
	@echo "install:           install NEMOBIN/gravidy"
	@echo "nbody6:            grab the code from URL1=$(URL1)"
	@echo "nbody6++:          grab the code from URL2=$(URL2)"
	@echo "install_nbody6:    install in $(NEMOBIN)"
	@echo "install_nbody6++:  install in $(NEMOBIN)"

install:	install_nbody6 install_nbody6++ install_nbody6++.avx


nbody6:
	$(WGET) $(URL1)
	(mkdir nbody6 ; cd nbody6; tar zxf ../nbody6.tar.gz)


nbody6++: Nbody6ppGPU

Nbody6ppGPU:
	git clone $(URL2)

# just an example on a simple laptop
# note the order of building appears to matter !!!
install_nbody6:  nbody6
	(cd nbody6/GPU2 ;  make clean; make sse ; cp run/nbody7.sse $(NEMOBIN))
	(cd nbody6/Ncode ; make clean; make     ; cp nbody6         $(NEMOBIN))

# just an example on a simple laptop
# the avx option may not work in some XEON's, watch out for the outut of configure

install_nbody6++.avx:  Nbody6ppGPU
	(cd Nbody6ppGPU ; make clean ;\
	./configure --disable-gpu --disable-mpi ;\
	make ;\
	cp build/nbody6++.avx $(NEMOBIN) ;\
	cp build/nbody6++.avx $(NEMOBIN)/nbody6++)

# this one should work more generically
install_nbody6++:  Nbody6ppGPU
	(cd Nbody6ppGPU ; make clean ;\
	./configure --disable-gpu --disable-mpi --enable-simd=sse ;\
	make ;\
	cp build/nbody6++.sse $(NEMOBIN) ;\
	cp build/nbody6++.sse $(NEMOBIN)/nbody6++)

# example runs

NBODY = 1024
TCRIT = 1.0

p.snap:
	mkplummer p.snap $(NBODY) seed=123
	mknemo runbody6

# nbody6++.SSE mode (also the default)
run1:   p.snap
	time runbody6 p.snap run1 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 exe=nbody6++.sse
	snapprint run1/OUT3.snap key,x,y,z,vx,vy,vz times=$(TCRIT) | head -4

# nbody6 single processor
run2:   p.snap
	time runbody6 p.snap run2 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=1 exe=nbody6 
	snapprint run2/OUT3.snap key,x,y,z,vx,vy,vz times=1 | head -4

# nbody6 SSE multiprocessor
run3:   p.snap
	time runbody6 p.snap run3 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=1 exe=nbody7.sse
	snapprint run3/OUT3.snap key,x,y,z,vx,vy,vz times=1 | head -4

# nbody6++.AVX mode
run4:   p.snap
	time runbody6 p.snap run4 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 nbody6=1 exe=nbody6++.avx
	snapprint run3/OUT3.snap key,x,y,z,vx,vy,vz times=1 | head -4

#
run9:	p.snap
	time runbody4 p.snap run9 tcrit=$(TCRIT) deltat=0.1 dtadj=0.05 

#       for kicks, a softened gravity run
run0:
	time rungravidy p.snap run0 tcrit=$(TCRIT) deltat=0.1 eps=0.01

