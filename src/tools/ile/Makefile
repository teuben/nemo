# vix 13dec91

DESTPATH = /usr/local
CFLAGS = -g
CC = cc

all: ile

clean:
	-rm *.o
	-rm ile

install:
	install -c ile $(DESTPATH)/bin/ile
	install -c ile.1 $(DESTPATH)/man/man1/ile.1

ile: ile.o
	cc -o ile ile.o -ltermcap
