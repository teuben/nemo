#
include $(NEMOLIB)/makedefs

# need to override the overly complex versions that makedefs have ?
CFLAGS = -O
FFLAGS = -O

#   For cambridge, or starlib, sites
#SL=/star/lib
#LIBP = -L$(SL) `pgplot_link`

LIBP = -L$(PGPLOT_DIR) -lpgplot -L/usr/X11R6/lib -lX11
LIBX = -L/usr/X11R6/lib -lX11
LIBX = -L/usr/lib/x86_64-linux-gnu/ -lX11

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


#	these are benchmarks with the (X) windows display disabled
#       use inx0 and inpyth for the ones with a window
bench:
	time xtriple < inx0
	time ptriple < inpyth0

