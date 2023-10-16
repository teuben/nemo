
FC = gfortran
FFLAGS = -g -O2

#LIBS = -lpgplot -lpng -L/usr/X11R6/lib -lX11
LIBS = -lpgplot

BINS = epic5

DUMMY = EPIRV EPIDEN EPITEST

all:  $(BINS) $(DUMMY)

clean:
	rm -f $(BINS)

EPIRV:
	@echo Creating dummy file $@ for running epic5
	@touch $@

EPIDEN:
	@echo Creating dummy file $@ for running epic5
	@touch $@

EPITEST:
	@echo Creating dummy file $@ for running epic5
	@touch $@

install: $(BINS)
	cp $(BINS) $(NEMOBIN)

epic5: epic5.f
	$(FC) $(FFLAGS) -o epic5 epic5.f $(LIBS)

