      SUBROUTINE INPARS(nmax,n,eta,deltat,tcrit,eps2,reset)
C
C Input:
C   nmax     max number particles in workspace
C Output: (The input integration parameters)
C   n        number of bodies to be read
C   eta      accuracy parameter
C   deltat   output time interval
C   tcrit    length of integration
C   eps2     square of softening length
C
      INCLUDE 'fdefs.inc'
      INTEGER nmax,n,reset
CSED    The next statement can be modified with SED to toggle type
      DOUBLE PRECISION    eta,deltat,tcrit,eps2

      WRITE(6,*) 'Enter n,eta,deltat,tcrit,eps2,reset:'
      READ(5,*) n,eta,deltat,tcrit,eps2,reset
      IF (n.LE.nmax) RETURN

      WRITE(6,*) 'Recompile program with larger value of NMAX'
      WRITE(6,*) 'Current value of NMAX = ',nmax
      WRITE(6,*) 'whereas number of bodies to be read, N = ',n
      STOP

      END
C***********************************************************************
      SUBROUTINE INBODS (n, body, x0, x0dot)
C       input of initial conditions for integrations
      INCLUDE 'fdefs.inc'
      INCLUDE 'nmax.inc'
      INTEGER i,k,n
CSED    The next statement can be modified with SED to toggle type
      DOUBLE PRECISION body(*), x0(NDIM,*), x0dot(NDIM,*)

      DO i=1,n
         READ(5,*) body(i),(x0(k,i),k=1,NDIM),(x0dot(k,i),k=1,NDIM)
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE OUTBODS (body, x, a, step, i)
C   output particle stuff
      INCLUDE 'fdefs.inc'
      INCLUDE 'nmax.inc'
CSED    The next statement can be modified with SED to toggle type
      DOUBLE PRECISION    body, x(NDIM), a(NDIM), step
      INTEGER i,k

      WRITE(6,105) body,(x(k),k=1,NDIM),(a(k),k=1,NDIM),step, i
  105 FORMAT(f10.2,3x,3f10.2,3x,3f10.2,3x,f12.4,3x,i5)
      RETURN
      END
C***********************************************************************
      SUBROUTINE OUTENE (tnext, nsteps, e)
C
C   Output current timestep, number of steps performed so far and
C   total energy of system
C
      INCLUDE 'fdefs.inc'
CSED    The next statement can be modified with SED to toggle type
      DOUBLE PRECISION    tnext, e
      INTEGER nsteps

      WRITE(6,140) tnext,nsteps,e
  140 FORMAT(5x,'time =',f7.2,'  steps =',i6,' energy =',f10.4,/)
      RETURN
      END
