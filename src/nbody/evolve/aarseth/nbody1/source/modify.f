      SUBROUTINE MODIFY(KSTART)
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common1.h'
      EXTERNAL VERIFY
*
*
*       Read first, second or both lines (KSTART = 3, 4, 5).
      IF (KSTART.EQ.4) GO TO 10
*
*       Read new DELTAT, TNEXT, TCRIT, QE & KZ(J) (if > 0).
      READ (5,*)  DT, TN, TC, QE1, J, K
*
*       Set new parameters if corresponding input is non-zero.
      IF (DT.LE.0.0) THEN
          DT = DELTAT
      ELSE
          DT = DT*TCR0
      END IF
*
      TDUM = TIME
      IF (TN.LE.0.0) THEN
          TNEXT = MAX(TNEXT - DELTAT + DT,TDUM)
      ELSE
          TNEXT = MAX(TN*TCR0,TDUM)
      END IF
*
      DELTAT = DT
      IF (TC.GT.0.0) TCRIT = TC*TCR0
      IF (QE1.GT.0.0) QE = QE1
*
*       See whether any options should be changed.
      IF (J.GT.0) KZ(J) = K
*
      WRITE (6,5)  DELTAT, TCRIT, QE, J, K
    5 FORMAT (///,7X,'RESTART PARAMETERS:   DELTAT =',F7.3,
     &          '  TCRIT =',F7.1,'  QE =',1P,E9.1,'  KZ(',I2,') =',I2,/)
*
*       Read new ETA (if > 0 & KSTART = 4 or 5).
   10 IF (KSTART.GE.4) THEN
          READ (5,*)  ETA1
*
*       Check modification of integration parameter.
          IF (ETA1.GT.0.0) ETA = ETA1
*
          WRITE (6,15)  ETA
   15     FORMAT (/,7X,'RESTART PARAMETERS:   ETA =',F7.3,/)
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Save the new parameters on tape/disc in case a restart is needed.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
