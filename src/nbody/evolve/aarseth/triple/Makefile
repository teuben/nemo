#
include $(NEMOLIB)/makedefs

CFLAGS = -O
FFLAGS = -O
FC = f77

#SL=/star/lib
#LIBP = -L$(SL) `pgplot_link`

LIBP = -lpgplot -L/usr/X11/lib -lX11

LIBX = -L/usr/X11/lib -lX11

SOURCE = triple.f binary.f block.f data.f derqp.f difsy1.f impact.f \
         peri.f transf.f wait.f

OBJECTS = $(SOURCE:.f=.o)

XOBJS = moviex.o xmovie.o Xtc.o 
POBJS = moviep.o 

help:
	@echo triple: Three-body regularization program.
	@echo Within NEMO probably use:
	@echo '    make LIBP="$$YAPP_PGPLOT" ...'
	@echo to set the PGPLOT library, the LIBX may need a new value too:
	@echo "Current value:  LIBX=$(LIBX)"

install:
	@echo Nothing to install

all:	xtriple ptriple

xtriple: $(OBJECTS) $(XOBJS)
	$(FC) $(FFLAGS) $(OBJECTS) $(XOBJS) $(LIBX) -o xtriple

ptriple: $(OBJECTS) $(POBJS)
	$(FC) $(FFLAGS) $(OBJECTS) $(POBJS) $(LIBP) -o ptriple

clean:
	rm -f core *.o xtriple ptriple 

# don't remove the next few lines, they are essential for Sverre

printh:
	@- \rm -f TRIPLE.TEXT
	@cath $(SOURCE) > TRIPLE.TEXT
print:
	@- \rm -f TRIPLE.TEXT
	@cat $(SOURCE) > TRIPLE.TEXT


#	these are benchmarks with the (X) windows display
bench:
	time xtriple < inx0
	time ptriple < inpyth0

