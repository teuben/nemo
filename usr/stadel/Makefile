#

URL1 = https://github.com/N-BodyShop/direct

URL2 = https://github.com/teuben/direct
BRANCH2 = teuben1

NBODY = 100
TIME = /usr/bin/time

direct:
	git clone -b $(BRANCH2) $(URL2)
	(cd direct; make)

install: direct
	(cd direct; cp direct $(NEMOBIN))

plummer.bin:
	mkplummer - $(NBODY) | snaptipsy - plummer.bin
	tipsysnap plummer.bin

bench: plummer.bin
	$(TIME) direct -o plummer < plummer.bin
	@echo plummer.acc and plummer.phi created

clean:
	rm -f plummer.???
