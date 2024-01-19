include $(NEMOLIB)/makedefs

URL1 = http://github.com/hannorein/rebound
URL2 = https://github.com/hannorein/reproducibility-paper
URL3 = https://github.com/dtamayo/reboundx

REBOUND = rebound
PROBLEM = simplest

help:
	@echo "install:         install rebound's lib and inc, test compile for PROBLEM=$(PROBLEM)"
	@echo "python:          install the python interface in your python"
	@echo "rebound:         grab the rebound code with as REBOUND=$(REBOUND) via github"
	@echo "pull:            ensure the code was updated to the latest"
	@echo "clean:           clean this directory"
	@echo "problems:        show available problems in rebound"

$(REBOUND):
	git clone $(URL1) $(REBOUND)

reproducibility-paper:
	git clone $(URL2)

pull: $(REBOUND)
	(cd $(REBOUND); git pull)


install: $(REBOUND)
	(cd $(REBOUND)/examples/$(PROBLEM); \
	make ; \
	cp librebound.so $(NEMOLIB); \
	cp ../../src/rebound.h $(NEMOINC))

python:	pip

pip: $(REBOUND)
	pip install -e $(REBOUND)

# -march=cpu-type      skylake-avx512
# cat /proc/cpuinfo | grep avx512

problems: $(REBOUND)
	@echo "The following example PROBLEM's are in rebound/examples:"
	@echo ""
	@(cd $(REBOUND)/examples; ls)
	@echo ""
	@echo "Use"
	@echo "   make install PROBLEM=..."
	@echo "to install one of the example problems"

clean:
	rm -rf $(REBOUND)

