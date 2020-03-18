#

WGET = wget
URL1 = http://www.h-its.org/wp-content/uploads/2014/11/galic.tar.gz
URL2 = https://wwwmpa.mpa-garching.mpg.de/~volker/galic/galic.tar.gz

GalIC:
	$(WGET) $(URL2)
	tar zxf galic.tar.gz
	cp GalIC/Template-Config.sh        GalIC/Config.sh
	cp GalIC/Template-Makefile.systype GalIC/Makefile.systype

NP = 4
M = `cd GalIC; echo Model_*.param`

models:
	for i in $(M); do \
	(cd GalIC; echo mpirun -np $(NP) GalIC $$i); done
