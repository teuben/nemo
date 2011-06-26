# ~pjt/orbit/potocode/Makefile: mantain source and executables for 
include $(NEMOLIB)/makedefs
# potcode

MBODY=16384
LM=-DMBODY=$(MBODY)

LDFLAGS =
L = $(NEMOLIB)/libnemo.a

# Objects exported to Nemo.

MAN1FILES = potcode.1
INCFILES = 
BINFILES = potcode
SRCFILES = defs.h potcode.c code_io.c orbstep.c dissipate.c diffuse.c
OBJFILES =        potcode.o code_io.o orbstep.o dissipate.o diffuse.o
SRCDIR = $(NEMO)/src/nbody/evolve/potcode
#
help:
	@echo No help for potcode yet

# Targets used to install potpole code in Nemo.

install: 

.install_man: $(MAN1FILES)
	cp $(MAN1FILES) $(NEMO)/man/man1
	@touch .install_man

.install_inc: 
	@echo Nothing 

.install_bin: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
	@touch .install_bin

.install_src: $(SRCFILES) $(INCFILES) Makefile
	@if [ ! -d $(SRCDIR) ]; \
	then \
		mkdir $(SRCDIR); \
		chmod 777 $(SRCDIR); \
	fi
	cp $? $(SRCDIR)
	@touch .install_src

# Targets used by Nemo to update include libraries and bin directories.

nemo_lib:
	@echo No nemo_lib here

nemo_bin: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)

bins: $(BINFILES)
	mv $(BINFILES) $(NEMOBIN)
#
clean:
	rm -f *.o *.BAK *.CKP potcode core

potcode: $(OBJFILES)
	$(CC) $(CFLAGS) -o potcode $(OBJFILES) $(NEMO_LIBS) $(FORLIBS)

potcode.o: potcode.c defs.h
	$(CC) $(LM) $(CFLAGS) -c potcode.c

orbstep.o: orbstep.c defs.h
	$(CC) $(LM) $(CFLAGS) -c orbstep.c

code_io.o: code_io.c defs.h
	$(CC) $(LM) $(CFLAGS) -c code_io.c

dissipate.o: dissipate.c defs.h
	$(CC) $(LM) $(CFLAGS) -c dissipate.c

diffuse.o: diffuse.c defs.h
	$(CC) $(LM) $(CFLAGS) -c diffuse.c


