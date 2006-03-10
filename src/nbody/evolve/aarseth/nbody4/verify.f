      SUBROUTINE VERIFY
*
*
*       Input validation.
*       -----------------
*
      INCLUDE 'common4.h'
*
*
*       Check for unreasonable input parameters (initial & restart).
      IF (N.GE.NMAX - 2) THEN
          WRITE (6,10)  N
   10     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5)
          STOP
      END IF
*
      IF (ETA.GT.0.08) THEN
          WRITE (6,20)  ETA
   20     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETA =',F6.2)
          STOP
      END IF
*
      IF (ETAU.GT.0.5.OR.GMIN.GT.0.0001.OR.GMAX.GT.0.10) THEN
          WRITE (6,30)  ETAU, GMIN, GMAX
   30     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAU =',F6.2,
     &                                 '  GMIN =',F11.7,'  GMAX =',F7.3)
          STOP
      END IF
*
*       Also check for zero or negative values.
      IF (N.LE.0.OR.ETA.LE.0.0) THEN
          WRITE (6,40)  N, ETA
   40     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5,
     &                                                   '  ETA =',F6.2)
          STOP
      END IF
*
      IF (ETAU.LE.0.0.OR.GMIN.LE.0.0.OR.GMAX.LE.0.0) THEN
          WRITE (6,30)  ETAU, GMIN, GMAX
          STOP
      END IF
*
      IF (DTADJ.LE.0.0.OR.DELTAT.LE.0.0.OR.QE.LE.0.0) THEN
          WRITE (6,50)  DTADJ, DELTAT, QE
   50     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   DTADJ =',F8.2,
     &                                '  DELTAT =',F6.2,'  QE =',1PE9.1)
          STOP
      END IF
*
      RETURN
*
      END
