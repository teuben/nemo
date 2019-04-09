include $(NEMOLIB)/makedefs

.KEEP_STATE:

#FFLAGS = -O3
#FC = f77


SOURCE = \
CGS.f init.f subroutines1.f legendre1.f diagnosticWORK.f plummer.f \
timestep.f arbel.f sfEraser.f f02aaf.f writepart.f


OBJECTS = $(SOURCE:.f=.o)

all:	CGS runCGS

install: all
	cp CGS $(NEMOBIN)/CGS.exe
	cp runCGS $(NEMOBIN)

CGS:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o CGS

bench:
	time ./CGS
	@echo "`tail -1 POS@EOS.DAT` `tail -1 VEL@EOS.DAT`"

bench2:
	rm -rf run_bench
	time runCGS run_bench
	tsf run_bench/snap.out | grep PhaseSpace | tail -1
	@echo Should look like:
	@echo "    double PhaseSpace[40000][2][3] -1.34270 0.860220 -2.62330 -0.197120 0.122590 0.0946220 -0.954800 -0.659210 -0.669670 0.380980 0.436340 0.286600 0.0282720 -0.283540"



clean:
	rm -f *.o CGS runCGS

