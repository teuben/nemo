
#

PY = amuse_convert.py amuse_integrate.py

help:
	@echo no help yet

amuse:
	git clone https://github.com/amusecode/amuse

git:	amuse
	(cd amuse; git pull)

install:
	cp $(PY) $(NEMOBIN)

profile:
	ipython profile create amuse
	@echo cp 95-amuse.py ~/.ipython/profile_amuse/startup/


ROW0 = "0.0078125 0.0 -0.256557735045 0.0863119421344 -0.0484937398944 -0.609540735858 -0.854122523728 0.736476012616"
ROW1 = "0.0078125 0.0 -0.147745056987 0.0102334055017 -0.107027598348 0.157821535395 -0.419219639057 -0.245930249344"

test1:
	./amuse_integrate.py 
	./amuse_convert.py -i plummer.amuse -I amuse -o plummer_0.txt -t 0
	./amuse_convert.py -i plummer.amuse -I amuse -o plummer_1.txt -t 10
	@echo "Phase space of first particle of first and last snapshot"
	@head -3 plummer_0.txt | tail -1
	@head -3 plummer_1.txt | tail -1
	@echo "For reproducability you should see:"
	@echo "$(ROW0)"
	@echo "$(ROW1)"
