#
# if your system is running 4.2 or is a sun running
# a sun os before 4.0 you need to add -DBSD to CFLAGS
#

CC = cc

MCFLAGS = -g -fsingle -f68881 /usr/lib/f68881/libm.il -DSUN_4_1_1

HDRS= ../src/vogl.h

DOBJS = sun.o ps.o hpdxy.o tek.o

CFLAGS= -I../src $(MCFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $*.c 

all:	$(DOBJS)

$(DOBJS): $(HDRS)

clean:
	rm -f core $(DOBJS)

clobber:
	rm -f core *.o
	
