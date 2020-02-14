#

TIME = /usr/bin/time
NP = 4


install:
	@echo To install gadget2 make sure you have LAM/MPI installed
	@echo As well as:   mknemo fftw2 gsl hdf5
	@echo After this
	@echo "             mknemo gadget2"
	@echo should work. 

# benchmarks 
galaxy:
	lamboot ;\
	( cd Gadget;\
	  mkdir -p galaxy;\
	  $(TIME) mpirun -np $(NP) Gadget2.galaxy Gadget2/parameterfiles/galaxy.param > galaxy.log;\
	); \
	lamwipe
	ls Gadget/galaxy
