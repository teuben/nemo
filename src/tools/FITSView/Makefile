#	 makefile for FITSView - but slightly modified to live within NEMO
include $(NEMOLIB)/makedefs
#	movie is called fitsmovie here - since we have a movie already
#	use CC macro
OBJS= ql rgb rgb24 fitsmovie
CFLAGS= -O
CC=/bin/cc
GLIBS=-lsuntool -lsunwindow -lpixrect -lm
OLIBS=loadimage.o fits.o
#
help:
	@echo Makefile for OBJS=$(OBJS)

all	: $(OBJS)

ql	: ql.o loadimage.o fits.o rainbow.o writeTGA.o
	$(CC) $(CFLAGS) ql.o loadimage.o fits.o rainbow.o writeTGA.o -o ql $(GLIBS)

fitsmovie   : movie.o loadimage.o fits.o rainbow.o movie_help.h
	$(CC) $(CFLAGS) movie.o loadimage.o fits.o rainbow.o -o fitsmovie $(GLIBS)

rgb	: rgb.o loadimage.o fits.o
	$(CC) $(CFLAGS) rgb.o loadimage.o fits.o  -o rgb $(GLIBS)

rgb24	: rgb24.o loadimage.o fits.o writeTGA.o
	$(CC) $(CFLAGS) rgb24.o loadimage.o fits.o writeTGA.o -o rgb24 $(GLIBS)


comb_install :
	install ql ../../../bin
	install rgb  ../../../bin
	install movie ../../../bin
	install rgb24 ../../../bin


clean:
	rm -f *.o $(OBJS)

fits.o : fits.h fitskeywords.h

ql.o rgb.o rgb24.o loadimage.o writeTGA.o: ql.h
