# include $(NEMOLIB)/makedefs

FC=gfortran
CC=gcc

FFLAGS=-g
CFLAGS=-g

SRCC=func_va_arg.c 
SRCF=main_va_arg.F

OBJC := ${SRCC:.c=.o} #substitute .f90 with .o
OBJF := ${SRCF:.F=.o} #substitute .f90 with .o
#   Example Makefile for hello.c within NEMO environment
#CFLAGS = -g
#LDFLAGS = -g
BIN = big
SRC = big.f

#	catch default
help:
	@echo Testing fortran array sizes, example:
	@echo 'make clean big FC=g77 NA=100 NB=100 NC=1000 ; ./big ; size big'



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

%.o: %.F #wildcard rule, creation of *.o depends on *.F
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.c #wildcard rule, creation of *.o depends on *.c
	$(CC) $(CFLAGS) -o $@ -c $<

test_va_arg : $(OBJC) $(OBJF)
	$(FC) $(FFLAGS) -o $@ $(OBJC) $(OBJF)

test: test_va_arg
	@./test_va_arg && echo "va_arg est passed" || echo "va_arg test failed"

clean:
	rm -f *.o *.a core $(BIN)
#	non-standard worker targets
