#
#   See FPDS/doc/doc_tutorial_cpp_en.pdf for more information
#

URL1 = https://github.com/FDPS/FDPS
URL2 = https://github.com/FDPS/PIKG

.PHONY:  test1 test2 test3

FDPS:
	git clone $(URL1) 
PIKG:
	git clone $(URL2)

# April 17-19, 2023. Honoring Jun's 60th.
test1:  FDPS
	(cd FDPS/sample/c++/nbody; make clean; make; \
	cp nbody.out $(NEMOBIN); \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=1/xs)

test2:  FDPS
	(cd FDPS/sample/c++/nbody; make clean; make CFLAGS="-std=c++17 -O3 -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp"; \
	cp nbody.out $(NEMOBIN); \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=2/xs)

test3:  FDPS
	(cd FDPS/sample/c++/nbody; make clean; make CFLAGS="-std=c++17 -O3 -DPARTICLE_SIMULATOR_MPI_PARALLEL" CC="mpicxx"; \
	cp nbody.out $(NEMOBIN); \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=3/xs)

test4:  FDPS
	(cd FDPS/sample/c++/nbody; make clean; make "use_phantom_grape_x86=yes"; \
	cp nbody.out $(NEMOBIN); \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=1/xs)

test5:  FDPS
	(cd FDPS/sample/c++/nbody; make clean; make CFLAGS="-std=c++17 -O3 -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp" use_phantom_grape_x86="yes"; \
	cp nbody.out $(NEMOBIN); \
	/usr/bin/time ./nbody.out; cd result; rm -f jun60; cat 0*.dat | tabtos - jun60 time,nbody key,mass,pos,vel; \
	snapplot jun60 times=0,3,6,9 nxy=2,2 xrange=-6:6 yrange=-6:6 yapp=2/xs)


phantom:
	(cd FDPS/src/phantom_grape_x86/G5/newton/libpg5/; make clean; make; cp libpg5.a $(NEMOLIB))
