
FLOATING_POINT = -fsingle /usr/lib/libm.il
CFLAGS = -O4 $(FLOATING_POINT) -I../../src -I/usr/local/R5/include/X11

EX = xtballs xtlcube

all: $(EX)

xtballs: xtballs.o balls.o 
	cc $(CFLAGS) -o xtballs xtballs.o balls.o ../../hershey/src/libhershey.a ../../src/libvogl.a -lXt -lXaw -lXmu -lXext -lX11 -lm

xtlcube: xtlcube.o lcube.o
	cc $(CFLAGS) -o xtlcube xtlcube.o lcube.o ../../src/libvogl.a -lXt -lXaw -lXmu -lXext -lX11 -lm

clean:
	rm -f $(EX) *.o core

clobber:
	rm -f $(EX) *.o core
