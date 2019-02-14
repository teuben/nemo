
help:
	@echo "install:         install NEMOBIN/gravidy"
	@echo "gravidy:         grab the code"
	@echo "ubuntu:          install BOOST (ubuntu)"
	@echo "gravidy-view:    install the 3D viewer"

gravidy:
	git clone http://github.com/cmaureir/gravidy.git

clean:
	rm -rf gravidy

install: gravidy
	(cd gravidy/src ; make pn; cp gravidy-cpu $(NEMOBIN)/gravidy)

gravidy-view:
	git clone https://gitlab.com/cmaureir/gravidy-view

ubuntu:
	sudo apt install libboost-all-dev


NEMO_BFLAGS = "BOOSTFLAGS =-I$(NEMO)/opt/include -I$(NEMO)/opt/lib -lboost_program_options"
install_gpu: gravidy
	(cd gravidy/src ; make gpu $(NEMO_BFLAGS); cp gravidy-cpu $(NEMOBIN)/gravidy)	
