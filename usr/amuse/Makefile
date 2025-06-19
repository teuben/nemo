#
TIME = /usr/bin/time
PY   = amuse_convert.py amuse_integrate.py
URL1 = https://github.com/amusecode/amuse

.PHONY : help

## help:       This Help
help : Makefile
	@sed -n 's/^##//p' $<

## amuse:      clone the code here (deprecated, use 'mknemo amuse')
amuse:
	git clone $(URL1)

git:	amuse
	(cd amuse; git pull)

## install:    scripts to NEMOBIN
install:
	cp $(PY) $(NEMOBIN)

## profile:    help making an ipython profile for AMUSE (ipython --profile amuse)
profile:
	ipython profile create amuse
	@echo cp 95-amuse.py ~/.ipython/profile_amuse/startup/


ROW0 = "0.0078125 0.0 -0.256557735045 0.0863119421344 -0.0484937398944 -0.609540735858 -0.854122523728 0.736476012616"
ROW1 = "0.0078125 0.0 -0.147745056987 0.0102334055017 -0.107027598348 0.157821535395 -0.419219639057 -0.245930249344"

NBODY = 2048

## test1:      basic test if amuse installed
test1:
	./amuse_integrate.py -o plummer.amuse
	./amuse_convert.py   -i plummer.amuse -o plummer_0.txt -t  0
	./amuse_convert.py   -i plummer.amuse -o plummer_1.txt -t 10
	@echo "Phase space of first particle of first and last snapshot"
	@head -3 plummer_0.txt | tail -1
	@head -3 plummer_1.txt | tail -1
	@echo "For reproducability you should see:"
	@echo "$(ROW0)"
	@echo "$(ROW1)"

## bench1a:    bench1 AMUSE make plummer
bench1a:
	$(TIME) ./amuse_convert.py -n 10000000 -s 123 
	$(TIME) ./amuse_convert.py -n 10000000 -s 123 -o big.amuse

## bench1n:    bench1 NEMO mkplummer
bench1n:
	@rm -f big.nemo
	$(TIME) mkplummer .           10000000
	$(TIME) mkplummer big.nemo    10000000

## bench2a:    bench2 AMUSE bhtree
bench2a:
	$(TIME) ./amuse_integrate.py -o sim.amuse -n $(NBODY) -s 123 -T 10
	@./amuse_convert.py -i sim.amuse -o plummer0.txt -t 0

## bench2n:    bench2 NEMO hackcode1
bench2n:
	@echo "Make sure bench2a was run before bench2n"
	@rm -f plummer0.dat plummer1.dat
	tabcols plummer0.txt | tabtos - plummer0.dat block1=mass,skip,vx,vy,vz,x,y,z nbody=$(NBODY)
	$(TIME) hackcode1 plummer0.dat plummer1.dat tol=0.75 tstop=10 freq=64 eps=0.025 > /dev/null
	@snapprint plummer1.dat times=10 header=nbody,time newline=t format=%16.13f >plummer1.tab


# pick bhtree or BHTC (see usr/kawai/nbody)
BHTREE = bhtree

## bench2b:    bench2 NEMO/bhtree
bench2b:
	@mkplummer - $(NBODY) | snapprint - header=nbody,time newline=t format=%16.13f >p2k.tab
	@rm -f p2k.out p2k.dat
	$(TIME) $(BHTREE) -i p2k.tab -o p2k.out -D 1 -T 10 -d 0.015625      >p2k.log 2>p2k.err
	tail p2k.err
	@atos p2k.out p2k.dat mass,phase,phi 
