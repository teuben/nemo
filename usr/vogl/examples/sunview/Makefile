#
# examples makefile 
#


#
# Where to find library
LIB = ../../hershey/src/libhershey.a ../../src/libvogl.a
LIBS = -lX11 -lsuntool -lsunwindow -lpixrect -lm
#
MCFLAGS = -g -fsingle /usr/lib/libm.il
CFLAGS = -I../../src $(MCFLAGS)

EX = svlcube svballs

all: $(EX)

svlcube: svlcube.o lcube.o $(LIB)
	$(CC) $(CFLAGS) -o svlcube svlcube.o lcube.o $(LIB) $(LIBS)

svballs: svballs.o balls.o $(LIB)
	$(CC) $(CFLAGS) -o svballs svballs.o balls.o $(LIB) $(LIBS)

clean:
	rm -f *.o core

clobber:
	rm -f $(CEXAMPS) *.o core $(EX)
