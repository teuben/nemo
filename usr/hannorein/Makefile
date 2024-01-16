include $(NEMOLIB)/makedefs

REBOUND = rebound
PROBLEM = simplest

help:
	@echo "install:         install NEMOBIN/rebound for PROBLEM=$(PROBLEM)"
	@echo "rebound:         grab the rebound code with as REBOUND=$(REBOUND) via github"
	@echo "pull:            ensure the code was updated to the latest"
	@echo "clean:           clean this directory"
	@echo "problems:        show available problems in rebound"

$(REBOUND):
	git clone http://github.com/hannorein/rebound   $(REBOUND)

pull: $(REBOUND)
	(cd $(REBOUND); git pull)


install: $(REBOUND)
	(cd $(REBOUND)/examples/$(PROBLEM); make ; cp librebound.so $(NEMOLIB); \
		cp rebound $(NEMOBIN); cp rebound $(NEMOBIN)/rebound_$(PROBLEM))
	@echo "Note: in NEMO you can also execute this as rebound_$(PROBLEM)"

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

