      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'commonp.h'
      INTEGER  nxtlst(NMAX)
*
*
*       Find all particles in next block (TNEXT = TMIN).
    1 TMIN = 1.0D+10
      DO 10 I = 1,N
          TMIN = MIN(TNEXT(I),TMIN)
   10 CONTINUE
*
      NXTLEN = 0
      DO 20 I = 1,N
          IF (TNEXT(I).EQ.TMIN) THEN
              NXTLEN = NXTLEN + 1
              NXTLST(NXTLEN) = I
          END IF
   20 CONTINUE
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
      LI = 0
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,30) I, NAME(I), NSTEPS, TIME, STEP(I), TIME/STEP(I)
*  30     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check next output time.
      IF (TIME.GT.TPRINT) THEN
          TIME = TPREV
          CALL OUTPUT
          GO TO 1
      END IF
*
*       Predict all coordinates & velocities to order FDOT.
      DO 40 J = 1,N
          S = TIME - T0(J)
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
          XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
          XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
          XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   40 CONTINUE
*
*       Save new time and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
*
*       Advance the pointer (<= NXTLEN) and select next particle index.
   50 LI = LI + 1
      IF (LI.GT.NXTLEN) THEN
          DO 35 L=1,NXTLEN
              I = NXTLST(L)
              DO 32 K = 1,3
                  X(K,I) = X0(K,I)
                  XDOT(K,I) = X0DOT(K,I)
   32         CONTINUE
   35     CONTINUE
          GO TO 1
      END IF
      I = NXTLST(LI)
      TIME = T0(I) + STEP(I)
*
*       Perform the integration step for body #I.
      CALL NBINT(I)
*
*       Increase counters and check timer & optional COMMON save.
      NSTEPS = NSTEPS + 1
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.100000) GO TO 50
*
      NTIMER = 0
*     IF (NSTEPS.GE.100000) THEN
*         NSTEPS = 0
*         IF (KZ(1).GT.1) CALL MYDUMP(1,1)
*     END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 50
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME, TCOMP, CPUTOT/60.0, ERRTOT, DETOT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',1P,E9.1,
     &                 '  TCOMP =',0P,F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6)
      END IF
*
      STOP
*
  100 RETURN
*
      END
