FFLAGS = -O
FORTRAN = f77

SOURCE = \
firstn.f data.f energy.f force.f \
input.f intgrt.f output.f ran2.f scale.f

OBJECTS = $(SOURCE:.f=.o)

firstn:	$(OBJECTS)
	$(FORTRAN) $(FFLAGS) $(OBJECTS) -o firstn
