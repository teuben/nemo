# @(#)Makefile	1.1 (TRW) 1/14/86
CFLAGS = -O

.c.o:
	${CC} ${CFLAGS} -c $*.c
	-ld -x -r $*.o
	mv a.out $*.o

all: libglob.a

libglob.a: glob.o
	ar cr libglob.a glob.o
	ranlib libglob.a
	chmod 644 libglob.a

install: all

clean:
	rm -f *.o *.a

