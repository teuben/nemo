include $(NEMOLIB)/makedefs

.KEEP_STATE:

#FFLAGS = -O3
#FC = f77


SOURCE = \
CGS.f init.f subroutines1.f legendre1.f diagnosticWORK.f plummer.f \
timestep.f arbel.f sfEraser.f f02aaf.f writepart.f


OBJECTS = $(SOURCE:.f=.o)

all:	CGS

install: all
	cp CGS $(NEMOBIN)/CGS.exe

CGS:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o CGS

bench:
	time ./CGS
	@echo "`tail -1 POS@EOS.DAT` `tail -1 VEL@EOS.DAT`"

clean:
	rm -f *.o CGS

