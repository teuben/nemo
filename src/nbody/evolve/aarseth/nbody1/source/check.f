      SUBROUTINE CHECK(DE)
*
*
*       Error check and restart.
*       ------------------------
*
      INCLUDE 'common1.h'
*
*
*       Set absolute value of error & time-step modification factor.
      DE = ABS(DE)
      ETACOR = 1.0
*
*       Check restart for large errors (two attempts permitted).
      IF (DE.LT.5.0*QE) GO TO 30
*
      IF (KZ(2).LE.1.OR.NDUMP.GE.2) THEN
          WRITE (6,10)
   10     FORMAT (/,9X,'CALCULATIONS HALTED * * *')
*       Increase NDUMP to prevent 3rd restart (safety check in routine MAIN).
          NDUMP = NDUMP + 1
          IF (KZ(1).GT.0.AND.KZ(2).GE.1) CALL MYDUMP(1,1)
          STOP
      END IF
*
*       Repeat the previous interval with reduced time-step parameters.
      TCOMP = CPU
      NTEMP = NDUMP
      CALL MYDUMP(0,2)
      CPU = TCOMP
      ETACOR = 0.5
      ETA = ETACOR*ETA
      NDUMP = NTEMP + 1
*       Control variable NDUMP used to prevent a third restart.
      WRITE (6,20)  TIME, ETA
   20 FORMAT (/,9X,'RESTART * * *   TIME =',F8.2,'  ETA =',F7.3)
      GO TO 50
*
*       Reset counter and check optional modification of accuracy parameters.
   30 NDUMP = 0
      IF (KZ(11).EQ.0) GO TO 50
*
      IF (DE.GT.QE) THEN
*       Continue calculation but reduce the time-step parameters.
          ETACOR = SQRT(QE/DE)
          ETA = ETACOR*ETA
          WRITE (6,40)  ETA
   40     FORMAT (8X,'ETA =',F7.3)
      ELSE IF (DE.LT.0.2*QE) THEN
*       Continue calculation with increased time-step parameters (< ETA0).
          IF (TIME.GT.0.0D0) THEN
              FAC = 1.2
              ETACOR = MIN(FAC,ETA0/ETA)
              ETA = ETACOR*ETA
              IF (ETACOR.GT.1.01) WRITE (6,40)  ETA
          END IF
      END IF
*
*       See whether the time-steps should be changed.
   50 IF (ETACOR.NE.1.0) THEN
          ETACOR = SQRT(ETACOR)
          DO 60 I = 1,N
              STEP(I) = ETACOR*STEP(I)
*       Check that time will not decrease.
              IF (T0(I) + STEP(I).LT.TIME) STEP(I) = TIME - T0(I)
   60     CONTINUE
      END IF
*
      RETURN
*
      END
