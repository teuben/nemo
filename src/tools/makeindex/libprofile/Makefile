# @(#)Makefile	1.1 (TRW) 1/14/86
CFLAGS = -O
OBJECTS = boolean.o free.o has.o read.o space.o write.o

.c.o:
	${CC} ${CFLAGS} -c $*.c
	-ld -r -x $*.o
	mv a.out $*.o

all: libprofile.a

libprofile.a: ${OBJECTS}
	ar cr libprofile.a ${OBJECTS}
	ranlib libprofile.a
	chmod 644 libprofile.a

install: all

clean:
	rm -f *.o *.a
