#   Makefile for the menu program, mirtool
include $(NEMOLIB)/makedefs
#
#	Intended as a general tool for any project with programs
#	of the execution nature 
#		'program [key=value ...]'
#
#	Last mod:   15-jul-90   all files in one exectool.c	PJT
#
#  change install directory IDIR while calling make:
#  e.g. 'make src IDIR=$MIR/src/menu/exectool'
IDIR = /tmp

#
CC=/bin/cc
CFLAGS	= -O -Bdynamic
SUNTOOLS= -lsuntool -lsunwindow -lpixrect

#	for now it's a kludge, one can make mirtool and nemotool
#	at the same time

# OLDSRCS	= main.c search.c lexec.c hsearch.c load.c
SRC = mirtool.c 

ALL	= $(SRC) Makefile mirtool.icon README*

help:
	@echo Makefile for mirtool and nemotool

mirtool:	mirtool.c
	$(CC) $(CFLAGS) -DMIR -o mirtool mirtool.c $(SUNTOOLS)

nemotool:	mirtool.c
	$(CC) $(CFLAGS) -DNEMO -o nemotool mirtool.c $(SUNTOOLS)

exectool:
	@echo "No such program, use mirtool or nemotool"
	@echo "Exectool is the generic name"

clean:	cleanup

cleanup:
	rm -f *.o nemotool mirtool *.def

#	Below this Old and Private stuff 
#
lint:
		lint $(SRC) | egrep -v "inconsist"

main_m.o:	main.c
	$(CC) $(CFLAGS) -c main.c
	mv main.o main_m.o

search_m.o:	search.c
	$(CC) $(CFLAGS) -c search.c
	mv search.o search_m.o

lexec_m.o:	lexec.c
	$(CC) $(CFLAGS) -c lexec.c
	mv lexec.o lexec_m.o

hsearch_m.o:	hsearch.c
	$(CC) $(CFLAGS) -c hsearch.c
	mv hsearch.o hsearch_m.o

load_m.o:	load.c
	$(CC) $(CFLAGS) -c load.c
	mv load.o load_m.o

main_n.o:	main.c
	$(CC) $(CFLAGS) -c main.c
	mv main.o main_n.o

search_n.o:	search.c
	$(CC) $(CFLAGS) -c search.c
	mv search.o search_n.o

lexec_n.o:	lexec.c
	$(CC) $(CFLAGS) -c lexec.c
	mv lexec.o lexec_n.o

hsearch_n.o:	hsearch.c
	$(CC) $(CFLAGS) -c hsearch.c
	mv hsearch.o hsearch_n.o

load_n.o:	load.c
	$(CC) $(CFLAGS) -c load.c
	mv load.o load_n.o


#	special install script, which copies the source to an
#	install directory

src:	$(ALL)
	cp $? $(IDIR)

install:  .install_src

submit:  .submit

.submit:  $(SRC)
	submit $? $(MIR)/src/tools
	@touch .submit

.install_src:  $(SRC) Makefile README.PJT
	cp $? $(NEMO)/src/tools/miriad/mirtool
	@touch .install_src

