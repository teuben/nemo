#  NEMO/src/nbody/evolve/flowcode

include $(NEMOLIB)/makedefs

MBODY=16384

LOCAL_INC =  -DMBODY=$(MBODY)

#CFLAGS = -g
#LDFLAGS =
L = $(NEMOLIB)/libnemo.a

# Objects exported to Nemo.

MAN1FILES = flowcode.1
INCFILES = 
BINFILES = flowcode
SRCFILES = defs.h flowcode.c code_io.c orbstep.c dissipate.c diffuse.c vxy.c
OBJFILES =        flowcode.o code_io.o orbstep.o dissipate.o diffuse.o
SRCDIR = $(NEMO)/src/nbody/evolve/flowcode
OBJ = vrt.so  vrtd.so  vxy.so

TARFILES = $(MAN1FILES) $(SRCFILES) Makefile

#
# code_io.c    diffuse.c    flowcode.c   vxy.c
# defs.h       dissipate.c  orbstep.c
#

help:
	@echo No help for flowcode yet

# Targets used to install flowcode code in Nemo.

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
	rm -f *.o *.BAK *.CKP *~ core flowcode 

tar:
	tar cvf flowcode.tar $(TARFILES)

# 	$(CC) $(CFLAGS) $(LOCAL_INC) -o $* $*.o $(NEMO_LIBS) $(LOCAL_LIB) $(FORLIBS) $(EL)

flowcode: $(OBJFILES)
	$(CC) $(CFLAGS) -o flowcode $(OBJFILES) $(NEMO_LIBS) $(FORLIBS)

flowcode.o: flowcode.c defs.h
	$(CC) $(CFLAGS) -c flowcode.c

orbstep.o: orbstep.c defs.h
	$(CC) $(CFLAGS) -c orbstep.c

code_io.o: code_io.c defs.h
	$(CC) $(CFLAGS) -c code_io.c

dissipate.o: dissipate.c defs.h
	$(CC) $(CFLAGS) -c dissipate.c

diffuse.o: diffuse.c defs.h
	$(CC) $(CFLAGS) -c diffuse.c


#




obj:  $(OBJ)

