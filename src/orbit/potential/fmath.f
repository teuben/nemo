c
c This routine should be compiled and linked with potential routines
c that need to load fortran compiled potentials and need math library
c routines. Usually means you should defined $FORLIBS for the C compiler
c 
c This routine is a maintenance nightmare.... one needs to make dummy
c calls to routines that you find any of the fortran potentials needs
c
	subroutine fmath
	double precision da,db
	real a
	a = 3
	a = a**1.5
	da = mod(da,db)
	end
