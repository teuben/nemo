# makefile for remotename

# locations on astro
BINDIR=/local/bin
MANDIR=/local/man

# locations on campus software
#BINDIR=/usr/local/bin
#MANDIR=/usr/share/local/man

CFLAGS=-O # -g -DDEBUG
LDFLAGS= #-g

remotename:	remotename.o
	${CC} -o $@ ${LDFLAGS} remotename.o

remotename.o:	remotename.c
	${CC} -c ${CFLAGS} remotename.c

debug:
	${MAKE} CFLAGS=-g LDFLAGS=-g

lint:
	lint remotename.c

install:	remotename FRC
	install -s remotename ${BINDIR}
	install -m 0644 remotename.1 ${MANDIR}/man1

uninstall:	FRC
	rm -f ${BINDIR}/remotename ${MANDIR}/man1/remotename.1
	

clean:
	rm -f remotename.o remotename OUT

FRC:
