      SUBROUTINE VERIFY
*
*
*       Input validation.
*       -----------------
*
      INCLUDE 'common2.h'
*
*
*       Check for unreasonable input parameters (initial & restart).
      IF (N.GT.NMAX.OR.NNBMAX.GE.LMAX) THEN
          WRITE (6,10)  N, NNBMAX
   10     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5,
     &                                                  '  NNBMAX =',I4)
          STOP
      END IF
*
      IF (ETAI.GT.0.04.OR.ETAR.GT.0.08) THEN
          WRITE (6,20)  ETAI, ETAR
   20     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAI =',F6.2,
     &                                                  '  ETAR =',F6.2)
          STOP
      END IF
*
*       Also check for zero or negative values.
      IF (N.LE.0.OR.NNBMAX.LE.0.OR.ETAI.LE.0.0.OR.ETAR.LE.0.0) THEN
          WRITE (6,40)  N, NNBMAX, ETAI, ETAR
   40     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5,
     &                  '  NNBMAX =',I4,'  ETAI =',F6.2,'  ETAR =',F6.2)
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
