      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common1.h'
*
*
*       Read & print the main input parameters.
      READ (5,*)  N, ETA, EPS, DELTAT, TCRIT
*
      WRITE (6,1)  N, ETA, EPS
    1 FORMAT (//,5X,'N =',I5,'  ETA =',F6.1,'  EPS =',F7.3)
*
*       Set square softening parameter and initialize times.
      EPS2 = EPS**2
      TIME = 0.0
      TNEXT = 0.0
      NSTEPI = 0
*
      RETURN
*
      END
