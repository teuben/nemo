c
c This routine should be compiled and linked with potential routines
c that need to load fortran compiled potentials and need math library
c routines. Usually means you should defined $FORLIBS for the C compiler
c 
	subroutine fmath
	real a
	a = 3
	a = a**1.5
	end
