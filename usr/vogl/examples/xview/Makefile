#
# examples makefile 
#


#
# Where to find library
LIB = ../../hershey/src/libhershey.a ../../src/libvogl.a 
#XV_LIBS = -L/usr/openwin/lib -lxview -lolgx
#XV_LIBS = -L/usr/local/R5/lib -lxview -lolgx
XV_LIBS = -lxview -lolgx
#
MCFLAGS = -g -fsingle 
CFLAGS = -I../../src $(XV_INCS) $(MCFLAGS) 

EX = xvlcube xvballs

all: $(EX)

xvlcube: xvlcube.o lcube.o $(LIB)
	$(CC) $(CFLAGS) -o xvlcube xvlcube.o lcube.o $(XV_LIBS) $(LIB) $(LIBS)

xvballs: xvballs.o balls.o $(LIB)
	$(CC) $(CFLAGS) -o xvballs xvballs.o balls.o $(XV_LIBS) $(LIB) $(LIBS)


clean:
	rm -f *.o core

clobber:
	rm -f *.o core $(EX)
