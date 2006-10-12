#
include $(NEMOLIB)/makedefs
#   Example Makefile for hello.c within NEMO environment
#CFLAGS = -g
#LDFLAGS = -g
BIN = big
SRC = big.f

#	catch default
help:
	@echo Testing fortran array sizes



NA = 10
NB = 10
NC = 10

big:	big.f
	@echo "C  -- automatically generated, do not edit"     > big.inc
	@echo "      INTEGER NA,NB,NC"                        >> big.inc
	@echo "      PARAMETER (NA=$(NA),NB=$(NB),NC=$(NC))"  >> big.inc
	@echo "      REAL array(NA,NB,NB)"                    >> big.inc
	@echo Created big.inc with NA,NB,NC
	$(FC) $(FFLAGS) -o big big.f

clean:
	rm -f *.o *.a core $(BIN)
#	non-standard worker targets
