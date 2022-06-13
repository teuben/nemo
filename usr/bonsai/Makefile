#

URL1 = https://github.com/treecode/Bonsai
URL2 = https://github.com/sciserver/Bonsai

Bonsai:
	git clone $(URL1)
	(cd Bonsai/runtime; \
	cmake -DUSE_CUB=0 -DUSE_DUST=0 -DUSE_MPI=0 -DUSE_MPIMT=0;\
	make)

bench:
	(ulimit -s unlimited; \
	cd Bonsai/runtime; \
	./bonsai2 --plummer 1000 --snapname bench --dev 0 --snapiter 1 -T 10)
