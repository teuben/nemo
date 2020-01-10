
PDIRS = read_eagle martini

all:
	@echo $(PDIRS) pull

pull: $(PDIRS)
	for i in $(PDIRS); do\
	(cd $$i; git pull); done

read_eagle:
	git clone https://github.com/jchelly/read_eagle

install: read_eagle
	(cd read_eagle; cmake -DCMAKE_INSTALL_PREFIX:PATH=$(NEMO) .; make; make install)

tools: martini
	(cd martini; python setup.py install)

martini:
	git clone https://github.com/kyleaoman/martini


#   installing should not use python setup.py
install2:
	pip3 install -e martini[eaglesource]

install3:
	pip3 install -e martini.git[eaglesource]
