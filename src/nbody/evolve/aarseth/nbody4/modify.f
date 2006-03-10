      SUBROUTINE MODIFY(KSTART)
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common4.h'
      EXTERNAL VERIFY
*
*
*       Read first, second or both lines (KSTART = 3, 4, 5).
      IF (KSTART.EQ.4) GO TO 10
*
*       Read new DTADJ, DELTAT, TADJ, TPRINT, TCRIT, QE & KZ(J) (if > 0).
      READ (5,*)  DTA, DT, TA, TN, TC, QE1, J, K
*
*       Set new parameters if corresponding input is non-zero.
      IF (DTA.LE.0.0) DTA = DTADJ
      IF (DT.LE.0.0) DT = DELTAT
      IF (TA.LE.0.0) TADJ = MAX(TADJ - DTADJ + DTA,TIME)
      IF (TN.LE.0.0) TPRINT = MAX(TPRINT - DELTAT + DT,TIME)
*
      DTADJ = DTA
      DELTAT = DT
      IF (TA.GT.0.0) TADJ = TIME + DTADJ
      IF (TN.GT.0.0) TPRINT = MAX(TIME,TN-TOFF) 
      IF (TC.GT.0.0) TCRIT = TC
      IF (QE1.GT.0.0) QE = QE1
*
*       See whether any options should be changed.
      IF (J.GT.0) KZ(J) = K
*
      WRITE (6,5)  DTADJ, DELTAT, TCRIT, QE, J, K
    5 FORMAT (///,7X,'RESTART PARAMETERS:   DTADJ =',F7.3,'  DELTAT =',
     &                            F7.3,'  TCRIT =',F7.1,'  QE =',1PE9.1,
     &                                            '  KZ(',I2,') =',I2,/)
*
*       Read new ETA, ETAU, DTMIN, RMIN, ncrit (if > 0 & KSTART = 4 or 5).
   10 IF (KSTART.GE.4) THEN
          READ (5,*)  ETA1, ETA2, DTM, RM, newncrit
*
*       Check modification of integration parameters.
          IF (ETA1.GT.0.0) ETA = ETA1
          IF (ETA2.GT.0.0) ETAU = ETA2
          IF (DTM.GT.0.0) THEN
              DTMIN = DTM
              SMIN = 2.0*DTM
          END IF
          IF (RM.GT.0.0) THEN
              RMIN = RM
              RMIN2 = RM**2
              RMIN22 = 4.0*RMIN2
          END IF
          if (newncrit.gt.0) ncrit = newncrit
*
          WRITE (6,15)  ETA, ETAU, DTMIN, RMIN, ncrit
   15     FORMAT (/,7X,'RESTART PARAMETERS:   ETA =',F7.3,'  ETAU =',
     &                        F6.2,'  DTMIN =',1PE9.1,'  RMIN =',E9.1,
     &         '  ncrit = ',i5/)
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Save the new parameters on tape/disc unit #1 just in case.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
