*       commonp.
*       -------
*
      PARAMETER (NMAX=150)
      IMPLICIT REAL*8  (A-H,O-Z)
*
      COMMON/NBODY/  X(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),F(3,NMAX),
     &               FDOT(3,NMAX),BODY(NMAX),XDOT(3,NMAX),D1(3,NMAX),
     &               D2(3,NMAX),D3(3,NMAX),STEP(NMAX),T0(NMAX),
     &               N,NMASS,NAME(NMAX),KZ(10),NSTEPS,NTIMER,NBLOCK,
     &               IDUM1,NDUM(4)
*
      COMMON/PARAMS/ CPU,CPU0,CPUTOT,ETA,DELTAT,TPRINT,TCRIT,QE,
     &               TWOPI,ONE3,ONE6,ONE9,ONE12,
     &               TIME,ZMASS,BE(3),CMR(4),CMRDOT(4),ZKIN,POT,
     &               ERRTOT,DETOT,TCR
*
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(40),TNEXT(NMAX)
