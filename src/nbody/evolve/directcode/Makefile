include $(NEMOLIB)/makedefs

L= $(NEMOLIB)/libnemo.a
LIBN = $(NEMO_LIBS)


BINFILES = directcode

SRCFILES = code.c code.h code_io.c defs.h grav.c util.c 

SRCDIR = $(NEMOPATH)/src/nbody/evolve/directcode

all:	$(BINFILES)

# Targets used to export files to Nemo.

install: .install_bin
#install:

.install_man: $(MANFILES)
	cp $? $(NEMOPATH)/man/man1
	@touch .install_man

.install_bin: $(BINFILES)
	cp $? $(NEMOBIN)

.install_src: $(SRCFILES) Makefile
	@if [ ! -d $(SRCDIR) ]; \
	then \
		mkdir $(SRCDIR); \
		chmod 777 $(SRCDIR) ; \
	fi
	cp $? $(SRCDIR)
	@touch .install_src

test:	$(BINFILES)

# Targets used by Nemo to mantain bin files.

nemo_lib:
	@echo no nemo_lib here

nemo_bin: $(BINFILES)
	mv $? $(NEMOBIN)
	rm -f *.o

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

nemo_src:
	-@for i in $(BINFILES); do \
	echo `pwd` $$i ; done 

clean:
	rm -f *.o *.a core $(BINFILES)

tidy:
	rm -f *.o $(BINFILES)
#
directcode: code.o code_io.o grav.o util.o
	$(CC) $(CFLAGS) -o directcode \
	   code.o code_io.o grav.o util.o $(LIBN) -lm

code.o: code.c defs.h code.h

code_io.o: code_io.c defs.h

util.o: util.c defs.h

load.o: load.c defs.h

grav.o: grav.c defs.h

util.o: util.c defs.h


bench:	bench1 bench2

bench1:
	time directcode nbody=1024 > /dev/null

bench2:
	time directcode nbody=10240 > /dev/null

