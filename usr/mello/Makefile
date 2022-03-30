
PDIRS = clustep galstep

all:
	@echo $(PDIRS) pull example_c example_g

pull: $(PDIRS)
	for i in $(PDIRS); do\
	(cd $$i; git pull); done

clustep:
	git clone https://github.com/elvismello/clustep

galstep:
	git clone https://github.com/elvismello/galstep

install: $(PDIRS)
	(cd clustep/clustep; make)
	@echo clustep is installed in place
	@echo galstep is already installed

example_c: clustep/clustep/snapshot_0000
	(cd clustep/clustep; python3 clustep.py -o snapshot_c)
	@ls -l clustep/clustep/snapshot_0000

example_g:
	(cd galstep/galstep; python3 galstep.py -i params_galaxy.ini -o snapshot_g)
	@ls -l galstep/galstep
