#	$(NEMOPATH)/src/tools/movie
# 	11-nov-90	NEMO 2.x new location
#
#  Since this directory is for SUN's only - a dynamic BFLAG will considerably
#  decrease size of executable. The $(CC) is often a script on SUN with
#  -Bstatic to aid loadobj() stuff.
BFLAGS = -Bdynamic
CC = cc 
CFLAGS = -Bdynamic -g -DNEMO
L = $(NEMOLIB)/libnemo.a
BIN=movie movie_sv
PIX=$(BFLAGS) -lpixrect
GRX=$(BFLAGS) -lsuntool -lsunwindow -lpixrect

#		Only to be executed by NEMO:

nemo_lib:
	@echo No nemo_lib here

nemo_bin:	$(BIN)
	mv $? $(NEMOBIN)

nemo_src:
	-@for i in $(BIN); do \
	echo `pwd` $$i ; done 
		 
#		remainder of Makefile for 'public' use

install:    .install_man .install_bin

.install_man:	movie.1
	cp movie.1 $(NEMOPATH)/man/man1
	@touch .install_man

.install_bin: $(BIN)
	cp $(BIN) $(NEMOBIN)

clean:
	rm -f *.o $(BIN) core

all:	$(BIN)

#	Private targets

movie: movie.c
	$(CC) $(CFLAGS) -o movie movie.c $(L) $(PIX) -lm

movie_sv: movie_sv.c
	$(CC) $(CFLAGS) -o movie_sv movie_sv.c $(L) $(GRX) -lm

