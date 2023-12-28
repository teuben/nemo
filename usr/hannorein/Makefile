include $(NEMOLIB)/makedefs

help:
	@echo "install:         install NEMOBIN/rebound"
	@echo "rebound:         grab the code"

rebound:
	git clone http://github.com/hannorein/rebound 

install:
	(cd rebound/examples/simplest; make ; ./rebound; \
		cp rebound $(NEMOLIB); cp librebound.so $(NEMOLIB))

# -march=cpu-type      skylake-avx512

clean:
	rm -rf rebound

