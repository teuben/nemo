.KEEP_STATE:
FFLAGS = -O3
#FFLAGS = -g -fnonstd
FC = f77
#FORTRAN = g77

SOURCE = \
hermit.f block.f cputim.f data.f energy.f \
fpoly1.f iblock.f input.f intgrt.f \
mydump.f nbint.f output.f ran2.f \
start.f stepi.f stepk.f steps.f \
tstep.f xvpred.f zero.f

OBJECTS = $(SOURCE:.f=.o)

all:	hermit

hermit:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o hermit

print:
	@- \rm -f HERMIT.TEXT
	@cat $(SOURCE) > HERMIT.TEXT
