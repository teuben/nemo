      SUBROUTINE NTINT(I)
*
*
*       Single star integration.
*       ------------------------
*
      INCLUDE 'common4.h'
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3)
*
*
*       Predict position and velocity up to end of STEP.
    1 J = I
      S = STEP(J)
      S1 = 1.5*S
      S2 = 2.0*S
      X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S + X0(1,J)
      X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S + X0(2,J)
      X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S + X0(3,J)
      XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
      XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
      XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*
*       Copy X & XDOT to scalars and initialize FIRR & FD.
      ITER = 0
      DO 5 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    5 CONTINUE
*
*       Evaluate the galactic force and first derivative.
      CALL XTRNLT(XI,XIDOT,FIRR,FD)
*
*       Include the corrector and set new T0, F, FDOT, D1, D2 & D3.
      DT = STEP(I)
      DTSQ = DT**2
      DT6 = 6.0D0/(DT*DTSQ)
      DT2 = 2.0D0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
      T0(I) = T0(I) + STEP(I)
*
      DO 10 K = 1,3
	  DF = 2.0D0*F(K,I) - FIRR(K)
	  FD6 = 6.0D0*FDOT(K,I)
	  SUM = FD6 + FD(K)
	  AT3 = 2.0D0*DF + DT*SUM
	  BT2 = -3.0D0*DF - DT*(SUM + FD6)
*
	  X0(K,I) = XI(K) + (0.6D0*AT3 + BT2)*DTSQ12
	  X0DOT(K,I) = XIDOT(K) + (0.75D0*AT3 + BT2)*DT13
*
*       Update the corrected values (OK for test particles).
          X(K,I) = X0(K,I)
          XDOT(K,I) = X0DOT(K,I)
*
          F(K,I) = 0.5D0*FIRR(K)
          FDOT(K,I) = ONE6*FD(K)
*
*       Form derivatives even though not needed for commensurate times.
	  D2(K,I) = (3.0D0*AT3 + BT2)*DT2
	  D3(K,I) = AT3*DT6
*       NOTE: These are real derivatives!
   10 CONTINUE
*
*       Specify new time-step by standard criterion.
      TTMP = TSTEP(FIRR,FD,D2(1,I),D3(1,I),ETA)
      DT0 = TTMP
      TTIME = T0(I)
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          TTMP = MIN(2.0*STEP(I),1.0D0)
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
      ELSE
          TTMP = STEP(I)
      END IF
*
*       Set new block step and update next time.
      STEP(I) = TTMP
      TNEXT(I) = STEP(I) + T0(I)
*
*       Increase step counter.
      NSTAIL = NSTAIL + 1
*
*       See whether to continue until end of large block-step.
      IF (TNEXT(I).LT.TIME) THEN
          ITER = ITER + 1
          IF (ITER.LT.10) GO TO 1
          WRITE (6,20)  I, TIME, STEP(I), DT0, FIRR
   20     FORMAT (' SMALL TIDAL STEP    I T DT DT0 F',I7,F8.3,1P,5E10.2)
      END IF
*
      RETURN
*
      END
