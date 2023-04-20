#
#   See FPDS/doc/doc_tutorial_cpp_en.pdf for more information
#

URL1 = https://github.com/FDPS/FDPS
URL2 = https://github.com/FDPS/PIKG

PDPS:
	git clone $(URL1)
PIKG:
	git clone $(URL2)

# April 17-19, 2023. Honering Jun's 60th.
test1:  PDPS
	(cd FDPS/sample/c++/nbody; make clean; make; \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=1/xs)

test2:  PDPS
	(cd FDPS/sample/c++/nbody; make clean; make CFLAGS="-std=c++17 -O3 -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp"; \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=2/xs)

test3:  PDPS
	(cd FDPS/sample/c++/nbody; make clean; make CFLAGS="-std=c++17 -O3 -DPARTICLE_SIMULATOR_MPI_PARALLEL" CC="mpicxx"; \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=3/xs)

