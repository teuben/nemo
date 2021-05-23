# an example of a rebound example in a local directory
# while keeping rebound in a standard NEMO location

export OPENGL=0

#RDIR = ../rebound/src
RDIR = $(NEMO)/usr/hannorein/rebound/src
include $(RDIR)/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(RDIR) -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(RDIR)
	@-rm -f librebound.so
	@ln -s $(RDIR)/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(RDIR) clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
