# ~josh/nbody/tools/Makefile: some simple N-body tools.
include $(NEMOLIB)/makedefs
#	22-nov-89	updated for HOSTTYPE

L = $(NEMOLIB)/libnemo.a

# Files exported to Nemo.

MAN1FILES = 

BINS  = pathlabel

SRCFILES = pathlabel.c

SRCDIR = $(NEMO)/src/josh/nbody/tools

BINFILES = pathlabel

# Targets used by me to mantain Nemo files.

#install: .install_man .install_bin .install_src .install_local

all:	$(BINFILES)

install:


.install_man: $(MAN1FILES)
	cp $? $(NEMOPATH)/man/man1
	@touch .install_man

.install_bin: $(BINS)
	cp $? $(NEMOBIN)
	@touch .install_bin

.install_src: $(SRCFILES) Makefile
	@if [ ! -d $(SRCDIR) ]; \
	then \
		mkdir $(SRCDIR); \
		chmod 777 $(SRCDIR) ; \
	fi
	cp $? $(SRCDIR)
	@touch .install_src

# Targets used by Nemo to mantain bin files.

nemo_lib:
	@echo no nemo_lib in `pwd`

nemo_bin: $(BINS)
	mv $? $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:
	-@for i in $(BINS); do \
	echo `pwd` $$i ; done 
		
		
clean:
	rm -f *.o *.a core $(BINS)

tidy:
	rm -f *.o $(BINS) $(LOCALBIN)
