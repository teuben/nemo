      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common1.h'
*
*
*       Find next body to be advanced and set new time.
    1 NNB = NLIST(1) + 1
      A1 = 1.0E+06
*
      DO 2 L = 2,NNB
          J = NLIST(L)
          A2 = T0(J) + STEP(J)
          IF (A2.LT.A1) THEN
              I = J
              A1 = A2
          END IF
    2 CONTINUE
*
      TIME = T0(I) + STEP(I)
      IF (TIME.LT.TLIST) GO TO 10
*
*       Form new time-step list and re-determine next body to be treated.
      NNB = 1
    3 TLIST = TLIST + DTLIST
*
      DO 4 J = 1,N
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              NNB = NNB + 1
              NLIST(NNB) = J
          END IF
    4 CONTINUE
*
      IF (NNB.EQ.1) GO TO 3
      NLIST(1) = NNB - 1
*
*       Modify time-step list interval to optimize the procedure.
      A2 = SQRT(FLOAT(N))/FLOAT(NLIST(1))
*
*       Include inertial factor to stabilize NLIST membership on SQRT(N).
      DTLIST = SQRT(A2)*DTLIST
      GO TO 1
*
*       Advance the integration step.
   10 CALL NBINT(I)
*
*       Obtain new time-step.
      CALL STEPI(I)
*
*       Check option for producing movie frames (maximum NFRAME).
      IF (KZ(7).GT.0) THEN
          IF (TIME.GT.TFRAME) THEN
              IFRAME = TIME/DELTAF
              IF (IFRAME.LT.NFRAME) THEN
                  CALL FRAME
              ELSE
                  TFRAME = TCRIT
              END IF
          END IF
      END IF
*
*       Check next output time.
      IF (TIME.GT.TNEXT) GO TO 50
*
*       Advance counter and check timer & optional COMMON save.
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.1000) GO TO 1
      NTIMER = 0
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,40)  TIME, TCOMP, CPUTOT, ERRTOT
   40     FORMAT (//,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F8.2,
     &                              '  CPUTOT =',F6.2,'  ERRTOT =',F9.5)
      END IF
*
      STOP
*
   50 RETURN
*
      END
