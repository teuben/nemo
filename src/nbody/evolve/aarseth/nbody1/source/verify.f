      SUBROUTINE VERIFY
*
*
*       Input validation.
*       -----------------
*
      INCLUDE 'common1.h'
*
*
*       Check for unreasonable input parameters (initial & restart).
      IF (N.GT.NMAX) THEN
          WRITE (6,10)  N
   10     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I4)
          STOP
      END IF
*
      IF (ETA.GT.0.04) THEN
          WRITE (6,20)  ETA
   20     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETA =',F6.2)
          STOP
      END IF
*
*       Also check for zero or negative values.
      IF (N.LE.0.OR.ETA.LE.0.0) THEN
          WRITE (6,40)  N, ETA
   40     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I4,
     &                                                   '  ETA =',F6.2)
          STOP
      END IF
*
      IF (DELTAT.LE.0.0.OR.QE.LE.0.0) THEN
          WRITE (6,50)  DELTAT, QE
   50     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   DELTAT =',F6.2,
     &                                                 '  QE =',1P,E9.1)
          STOP
      END IF
*
      RETURN
*
      END
