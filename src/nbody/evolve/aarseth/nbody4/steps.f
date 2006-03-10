***
      SUBROUTINE STEPS(I1,I2,KCASE)
*
*
*       Initialization of time-steps & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Set new steps and initialize prediction variables.
      DO 40 I = I1,I2
*
*       Obtain time-step using DT = 0.5*ETA*F/FDOT (SQRT(ETA/F) if FDOT = 0).
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
          FD2 = FDOT(1,I)**2 + FDOT(2,I)**2 + FDOT(3,I)**2
*       Include precaution for small velocities (i.e. DT = 0.5*ETA).
          IF (FD2.LT.FI2) FD2 = FI2
          DT = 0.5*ETA*SQRT(FI2/FD2)
          IF (N.LE.10) DT = 0.1*DT
*
*       Reduce step for special case of triple, quad, chain or merger.
          IF (IPHASE.GE.4) DT = 0.5*DT
*
*       Include lower limit.
          IF (DT.LT.6.0E-11.AND.IPHASE.NE.-2) THEN
              WRITE (6,5)  I, IPHASE, SQRT(FI2), SQRT(FD2), DT
    5         FORMAT (' WARNING!    STEP    I IPH F FD DT ',
     &                                      I6,I4,1P,3E10.2)
              DT = MAX(DT,DTK(39))
          END IF
*
*       Initialize the time and obtain discrete time-step if required.
          T0(I) = TIME
*
*       Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DT,DTN)
          IF (TIME.LE.0.0D0) THEN
              STEP(I) = DTN
          ELSE 
*       Reduce step by factor 2 until commensurate with current time.
              STEP(I) = DTN
              ITER = 0
   10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
                  STEP(I) = 0.5D0*STEP(I)
                  ITER = ITER + 1
                  IF (ITER.LT.16.OR.STEP(I).GE.DTK(40)) GO TO 10
                  STEP(I) = DTK(40)
                  WRITE (6,12) I, ITER, STEP(I)
   12             FORMAT (' DANGER!    STEP    I ITER STEP ',
     &                                         I6,I4,1P,E10.2)
*                 CALL gpfree
*                 STOP
              ELSE IF (ITER.GE.16) THEN
                  WRITE (6,15) I, ITER, (TIME+TOFF)/STEP(I), DT, STEP(I)
   15             FORMAT (' WARNING!    STEP    I ITER T/STEP DT STEP ',
     &                                  I6,I4,F16.4,1P,2E9.1)
              END IF
          END IF
*
*       Update array for new block times (#I copied to free location).
          TNEXT(I) = T0(I) + STEP(I)
*
*       Set prediction variables (X0DOT set by START, KSREG or KSTERM).
          DO 30 K = 1,3
              X0(K,I) = X(K,I)
              F(K,I) = 0.5D0*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
              D2(K,I) = 0.0
              D3(K,I) = 0.0
   30     CONTINUE
   40 CONTINUE
*
      RETURN
      END
***
