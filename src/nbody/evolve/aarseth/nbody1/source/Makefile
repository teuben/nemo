#
include $(NEMOLIB)/makedefs
#FFLAGS = -O
#FC = f77
#	some linux versions need $(NEMOLIB)/etime.o here
LIBS =

SOURCE = \
nbody1.f block.f bodies.f cmcorr.f \
check.f cputim.f data.f define.f energy.f \
escape.f fpoly1.f fpoly2.f frame.f \
input.f intgrt.f modify.f mydump.f \
nbint.f output.f ran2.f remove.f \
scale.f start.f stepi.f steps.f \
subsys.f verify.f xvpred.f zero.f

OBJECTS = $(SOURCE:.f=.o)

help:
	@echo nbody1:
	@echo LIBS=$(LIBS)

nbody1:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o nbody1 $(LIBS)

install: nbody1
	mv nbody1 $(NEMOBIN)

bins: nbody1
	mv nbody1 $(NEMOBIN)


print:
	@- \rm -f NBODY1.TEXT
	@cat $(SOURCE) > NBODY1.TEXT

clean:
	\rm -f core nbody1 *.o


nmax:   .dummy

PARAMS = params.h
NMAX = 100

.dummy:
	@echo "Creating new $(PARAMS) with NMAX=$(NMAX)"
	@echo "*       NBODY1 parameters." > $(PARAMS)
	@echo "*       ------------------" >> $(PARAMS)
	@echo "*" 				>> $(PARAMS)
	@echo "      PARAMETER  (NMAX=$(NMAX))" >> $(PARAMS)
	@echo "*" 				>> $(PARAMS)
	@echo "*" 				>> $(PARAMS)
	@echo "*       ----------------------------------------" >> $(PARAMS)
	@echo "*       NMAX    Maximum number of bodies." >> $(PARAMS)
	@echo "*       ----------------------------------------" >> $(PARAMS)
	@echo "*" 						>> $(PARAMS)

