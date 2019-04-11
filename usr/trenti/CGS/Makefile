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

DAT2 = fort.90

bench2:
	rm -rf run_bench
	time runCGS run_bench nemo=f
	@tail -1 run_bench/$(DAT2)
	@echo Should look like:
	@echo "      40000   0.33825E+00   0.22327E+00  -0.79689E+00  -0.13400E+01  -0.55406E+00   0.16339E+00  -0.12474E+01  -0.88467E+00   0.30056E+01  -0.42294E+01"


# recall CGS uses G=4.4971 so the plummer needs to be rescaled for virial equilibrium

bench3:
	rm -rf p10k run_p10k
	mkplummer - 10000 | snapscale - p10k vscale='sqrt(4.4971)'
	time runCGS run_p10k in=p10k
	snapmass run_p10k/snap.out - 1/10000.0 | snapmradii - | tabplot - 1 2:10 line=1,1


clean:
	rm -f *.o CGS runCGS

