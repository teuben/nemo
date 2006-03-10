      SUBROUTINE CHECK(DE)
*
*
*       Error check and restart.
*       ------------------------
*
      INCLUDE 'common4.h'
*
*
*       See whether output intervals should be increased (at main output).
      IF (KZ(32).GT.0.AND.TIME.GE.TPRINT - 20.0*DTMIN) THEN
*       Check current energy level (factor of 2) for possible increase.
          K = KZ(32)
          ECRIT = 0.25/2.0**K
          IF (ABS(E(3)).LT.ECRIT) THEN
*       Define dynamical crossing time in case energy is near zero.
              TDYN = 2.0*RSCALE/SQRT(2.0*ZKIN/ZMASS)
              IF (2.0*DTADJ.GT.TDYN.OR.TIME.EQ.0.0) GO TO 5
              DTADJ = 2.0*DTADJ
              DELTAT = 2.0*DELTAT
              QE = SQRT(2.0)*QE
              KZ(32) = KZ(32) + 1
              WRITE (6,1)  DTADJ, DELTAT, QE
    1         FORMAT (/,5X,'NEW INTERVALS:   DTADJ =',F6.2,
     &                     '  DELTAT =',F6.2,'  QE =',1P,E8.1)
          END IF
      END IF
*
*       Perform automatic error check & restart (option 2).
    5 DE = ABS(DE)
      ETACOR = 1.0
*
*       Check restart for large errors (two attempts permitted).
      IF (DE.LT.5.0*QE) GO TO 30
      ETACOR = MAX(0.5D0,SQRT(QE/DE))
*
*       Terminate run if no further restart is allowed.
      IF (KZ(2).LE.1.OR.NDUMP.GE.2) THEN
          WRITE (6,10)  DE, 5.0*QE
   10     FORMAT (/,5X,'CALC. HALTED  -  ENERGY TOLERANCE EXCEEDED ',
     &                 '  DE =',1P,E9.1,'  5*QE =',E8.1) 
*       Increase NDUMP to prevent 3rd restart (safety check in routine MAIN).
          NDUMP = 2
          IF (KZ(1).NE.0.AND.KZ(2).GE.1) CALL MYDUMP(1,1)
          CALL gpfree
          STOP
      END IF
*
*       Repeat the previous interval with reduced time-step parameters.
      TCOMP = CPU
      NTEMP = NDUMP
      CALL MYDUMP(0,2)
      CPU = TCOMP
      NDUMP = NTEMP + 1
*       Control variable NDUMP used to prevent a third restart.
      ETA = ETACOR*ETA
      IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
      DTMIN = SQRT(ETACOR)*DTMIN
      SMIN = SQRT(ETACOR)*SMIN
      WRITE (6,20)  TTOT, ETA, ETAU
   20 FORMAT (/,9X,'RESTART * * *   TIME =',F7.1,'  ETA =',F7.3,
     &                                           '  ETAU =',F7.3)
      CALL MYDUMP(1,2)
      GO TO 50
*
*       Reset counter and check optional modification of accuracy parameters.
   30 NDUMP = 0
      IF (KZ(17).EQ.0) GO TO 50
*
      IF (DE.GT.QE) THEN
*       Continue calculation but reduce the time-step parameters.
          ETACOR = SQRT(QE/DE)
          ETA = ETACOR*ETA
          IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
          DTMIN = SQRT(ETACOR)*DTMIN
          SMIN = SQRT(ETACOR)*SMIN
          IF (ETACOR.LT.0.95) WRITE (6,40)  ETA, ETAU
   40     FORMAT (8X,'ETA =',F7.3,'  ETAU =',F7.3)
      ELSE IF (DE.LT.0.2*QE) THEN
*       Increase the time-step parameters (up to initial value only).
          IF (TIME.GT.0.0D0) THEN
              ETACOR = MIN(1.2D0,ETA0/ETA)
              ETA = ETACOR*ETA
              IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
              DTMIN = SQRT(ETACOR)*DTMIN
              SMIN = SQRT(ETACOR)*SMIN
              IF (ETACOR.GT.1.05) WRITE (6,40)  ETA, ETAU
          END IF
      END IF
*
*       See whether the time-steps should be reduced.
   50 IF (ETACOR.LT.1.0.AND.KZ(2).GT.2) THEN
          ETACOR = SQRT(ETACOR)
          DO 60 I = IFIRST,NTOT
              IF (DMOD(T0(I),0.5D0*STEP(I)).EQ.0.0D0) THEN
                  IF (T0(I) + 0.5*STEP(I).GE.TIME) THEN
                      STEP(I) = 0.5*STEP(I)
                      TNEXT(I) = T0(I) + STEP(I)
                  END IF
              END IF
*       Safety precaution to avoid bunching of time-steps in NLIST.
          ISEND = -1
   60     CONTINUE
*
      END IF
*
      RETURN
*
      END
