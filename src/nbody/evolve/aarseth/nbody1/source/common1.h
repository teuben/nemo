*       common1.
*       --------
*
      INCLUDE 'params.h'
      REAL*8         X,X0,X0DOT,T0,TIME
*
      COMMON/NBODY/  X(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),T0(NMAX),
     &               F(3,NMAX),FDOT(3,NMAX),BODY(NMAX),XDOT(3,NMAX),
     &               D1(3,NMAX),D2(3,NMAX),D3(3,NMAX),
     &               STEP(NMAX),T1(NMAX),T2(NMAX),T3(NMAX)
*
      COMMON/NAMES/  N,NFIX,NPRINT,NDUMP,MODEL,NRUN,NTIMER,NSTEPI,
     &               NFRAME,KZ(15),NLIST(NMAX),NAME(NMAX)
*
      COMMON/PARAMS/ TIME,CPU,ETA,DELTAT,TCRIT,QE,EPS2,
     &               ONE3,ONE6,ONE9,ONE12,TCR0,ETA0,RTIDE,
     &               CPU0,CPUTOT,ERRTOT,POT,ZKIN,DELTAF,TFRAME,
     &               TNEXT,TLIST,DTLIST,ZMASS,RSCALE,TCR,BE(3),
     &               CMR(4),CMRDOT(4),RBAR,ZMBAR,VSTAR,TSTAR
