*       common1.
*       --------
*
      PARAMETER (NMAX=256)
      IMPLICIT REAL*8  (A-H,O-Z)
*
      COMMON/NBODY/  BODY(NMAX),X(3,NMAX),XDOT(3,NMAX),F(3,NMAX),
     &               F1(3,NMAX),a1(3,NMAX),a2(3,NMAX),b1(3,NMAX),
     &               b2(3,NMAX),d2(3,NMAX)
*
      COMMON/NAMES/  N,NSTEPI,ICL,JCL
*
      COMMON/PARAMS/ TIME,ETA,DELTAT,TNEXT,TCRIT,STEP,STEP0,STEP1,
     &               EPS2,ZMASS,ZKIN,POT,BE(3),CMR(4),CMRDOT(4)
