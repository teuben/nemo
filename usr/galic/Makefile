#

URL1 = http://www.h-its.org/wp-content/uploads/2014/11/galic.tar.gz
URL2 = https://wwwmpa.mpa-garching.mpg.de/~volker/galic/galic.tar.gz

GalIC:
	wget $(URL2)
	tar zxf galic.tar.gz

M = `cd GalIC; echo Model_*.param`
models:
	for i in $(M); do \
	(cd GalIC; echo mpirun -np 8 $$i); done
