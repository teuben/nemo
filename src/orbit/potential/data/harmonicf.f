C inipotential, potential
C
C------------------------------------------------------------------------------
C INIPOTENTIAL: initializes the potential.
C      input: npar, the number of input parameters
C             par[] an array of npar parameters	
C	      name  string with some other optional information
C      If npar=0 defaults are taken 
C	remember to initialize them in a DATA statement in INIPOTENTIAL()
C	a COMMON block will have to be used to share info with POTENTIAL()
C------------------------------------------------------------------------------
C
	SUBROUTINE inipotential (npar, par, name)
C subroutine variables
	INTEGER            npar
	DOUBLE PRECISION   par(*)
        CHARACTER          name*(*)
C local variables
	INTEGER i
C global data
	DOUBLE PRECISION   omega,h(3)
	COMMON  /potpar/omega,h
C 	* Next DATA statement is actually illegal in ANSI fortran...
	DATA	omega/0.0d0/
	DATA    h/1.0d0,1.0d0,1.0d0/
C
	IF (npar.gt.0) THEN
	   omega = par(1)
	ENDIF
	DO 10 i=2,npar
	   h(i-1) = par(i) * par(i)
 10	CONTINUE
	par(1) = omega
C
C  *** IN RUN MODE THE WRITE STATEMENTS MUST BE COMMENTED OUT ***
C  *** For sake of ease of the Fortran_to_C interface in NEMO ***
C	write (*,'(''INIPOTENTIAL Harmonic - fortran version'')')
C	write (*,'(''Parameters : Pattern Speed ='',G14.7,
C     -                          '' (should be forced 0.0) '')') omega
C	write (*,'(''   wx,wy,wz= '',3G14.7)') h(1),h(2),h(3)
C  *** IN RUN MODE THE WRITE STATEMENTS MUST BE COMMENTED OUT ***
C
	END
C
C
C------------------------------------------------------------------------------
C  POTENTIAL: the worker routine. Determines at any given point x,y,z) the
C      forces and potential. Although the naming of parameters suggests
C      cartesian coordinates, they need not necessarely be, as long as the
C      the equations of motion can be written in the same way.
C      Note that this routine is good for 2 and 3D
C------------------------------------------------------------------------------
C
	SUBROUTINE POTENTIAL (ndim,pos,acc,pot,time)
C subroutine variables
	INTEGER           ndim
	DOUBLE PRECISION  pos(*), acc(*), pot, time
C local variables
	INTEGER i
C global data
	DOUBLE PRECISION  omega,h(3)
	COMMON  /potpar/omega,h
C	
	pot = 0.0d0
	DO 10 i=1,ndim
	   pot = pot + h(i)* pos(i)*pos(i)
	   acc(i) = -h(i) * pos(i)
 10	CONTINUE
	pot = 0.5d0 * pot

	END
