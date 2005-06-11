include $(NEMOLIB)/makedefs

# adapted from 'treecode 1.4' for NEMO V3

#	choose from:  -DQUICKSCAN, -DUSEFREQ, -DBINARYIO -DUSE_NEMO_IO
TFLAGS = -DUSE_NEMO_IO

L= $(NEMO_LIBS)
LIBN = $(L)

BINFILES = treecode1 treecode1_q

#SRCFILES = code.c code.h code_io.c defs.h grav.c load.c util.c hackforce.c
SRCFILES = treecode.c  treegrav.c  treeio.c  treeload.c


SRCDIR = $(NEMO)/src/nbody/evolve/hackcode/treecode1

all:	$(BINFILES)

bins:	all
	mv $(BINFILES) $(NEMOBIN)

libs:
	@echo Nothing to libs in $(SRCDIR)	

install:
	@echo Nothing to install in $(SRCDIR)

clean:
	rm -f *.o *.a core $(BINFILES)

tidy:
	rm -f *.o $(BINFILES)

treecode1: treecode.o  treegrav.o  treeio.o  treeload.o
	$(CC) $(CFLAGS) $(TFLAGS) -o treecode1 \
	   treecode.o  treegrav.o  treeio.o  treeload.o \
	  $(LIBN) -lm

treecode1_q: treecode_q.o  treegrav_q.o  treeio_q.o  treeload_q.o
	$(CC) $(CFLAGS) -o treecode1_q \
	   treecode_q.o  treegrav_q.o  treeio_q.o  treeload_q.o \
	  $(LIBN) -lm

# --- 

treecode_q.o:	treecode.c
	$(CC) $(CFLAGS) $(TFLAGS) -DQUICKSCAN -c treecode.c -o treecode_q.o

treegrav_q.o:	treegrav.c
	$(CC) $(CFLAGS) $(TFLAGS) -DQUICKSCAN -c treegrav.c -o treegrav_q.o

treeio_q.o:	treeio.c
	$(CC) $(CFLAGS) $(TFLAGS) -DQUICKSCAN -c treeio.c -o treeio_q.o

treeload_q.o:	treeload.c
	$(CC) $(CFLAGS) $(TFLAGS) -DQUICKSCAN -c treeload.c -o treeload_q.o

# ---

treecode.o:	treecode.c
	$(CC) $(CFLAGS) $(TFLAGS) -c treecode.c -o treecode.o

treegrav.o:	treegrav.c
	$(CC) $(CFLAGS) $(TFLAGS) -c treegrav.c -o treegrav.o

treeio.o:	treeio.c
	$(CC) $(CFLAGS) $(TFLAGS) -c treeio.c -o treeio.o

treeload.o:	treeload.c
	$(CC) $(CFLAGS) $(TFLAGS) -c treeload.c -o treeload.o

# ---
