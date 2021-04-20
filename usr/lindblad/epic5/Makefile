
FC = gfortran
FFLAGS = -g -O2

LIBS = -lpgplot

BINS = epic5

all:  $(BINS)

clean:
	rm -f $(BINS)

epic5:
	$(FC) $(FFLAGS) -o epic5 epic5.f $(LIBS)
