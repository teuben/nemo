
FC = gfortran
FFLAGS = -g -O2

#LIBS = -lpgplot -lpng -L/usr/X11R6/lib -lX11
LIBS = -lpgplot

BINS = epic5

all:  $(BINS)

clean:
	rm -f $(BINS)

epic5: epic5.f
	$(FC) $(FFLAGS) -o epic5 epic5.f $(LIBS)
