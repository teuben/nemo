# $Header$
LDFLAGS=-s
CFLAGS=-O -I../libglob -I../libprofile -I../include.local

all: index

index: indexlex.c index.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o index index.o ../getopt.o ../shift.o \
	../libprofile/libprofile.a ../libglob/libglob.a -ll

indexlex.c: indexlex.l

index.o: indexlex.c

install: all

clean:
	-rm *.o lex.yy.c y.tab.c indexlex.c index
