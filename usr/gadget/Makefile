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
	( cd Gadget;\
	  mkdir -p galaxy;\
	  $(TIME) mpirun -np $(NP) Gadget2.galaxy Gadget2/parameterfiles/galaxy.param > galaxy.log;\
	)
	ls Gadget/galaxy
	@echo Save the last snapshot in galaxy_snap.last
	gadgetsnap Gadget/galaxy/snapshot_007 Gadget/galaxy/galaxy_snap.last

# compare the results with another integrator in NEMO
galaxy_bench:
	rm -f Gadget/galaxy/snap0.in Gadget/galaxy/snap0.out
	gadgetsnap Gadget/galaxy/snapshot_000 - | snapscale - Gadget/galaxy/snap0.in mscale=43007.1
	gyrfalcON Gadget/galaxy/snap0.in Gadget/galaxy/snap0.out 3 0.5 kmax=8 eps=0.1
