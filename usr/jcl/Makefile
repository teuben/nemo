#

PDIRS = unsio uns_projects glnemo2

all:
	@echo $(PDIRS) pull

unsio:
	git clone https://gitlab.lam.fr/jclamber/unsio

uns_projects:
	git clone https://gitlab.lam.fr/jclamber/uns_projects

glnemo2:
	git clone https://gitlab.lam.fr/jclamber/glnemo2.git


pull: $(PDIRS)
	for i in $(PDIRS); do\
	(cd $$i; git pull); done
