FFLAGS = -O2
FC = f77

SOURCE = \
firstn.f data.f energy.f force.f \
input.f intgrt.f output.f ran2.f scale.f

OBJECTS = $(SOURCE:.f=.o)

all:	firstn

firstn:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o firstn

clean:
	rm -f *.o firstn core

N = 16 32 64 128 256

bench: firstn
	-@for i in $(N); do\
        (time firstn < input$$i > bench$$i.log); done
 

