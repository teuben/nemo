      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'commonp.h'
*
*
*       Options.
*       --------------------------------------------------
*       1   Common dump on unit #1 on'touch STOP'.
*       2   Common dump on unit #2 at TCRIT.
*       3   
*       4 
*       5   Alternative time-step expression. 
*       6
*       --------------------------------------------------
*
*       Read & print the main input parameters.
      READ (5,*)  N, NRAND, ETA, DELTAT, TCRIT
      READ (5,*)  (KZ(J),J=1,10)
      NMASS = N
*
      WRITE (6,5)  N, ETA
    5 FORMAT (//,5X,'N =',I3,'  ETA =',F8.5)
      WRITE (6,10)  (KZ(J),J=1,10)
   10 FORMAT (/,5X,'OPTIONS:   ',10I4,/)
*
*       Convert output times and decay time from years to scaled units.
*     DELTAT = TWOPI*DELTAT
*     TCRIT = TWOPI*TCRIT
*
*       Set random number skip for routine DATA.
      IDUM1 = NRAND
*
      RETURN
*
      END
