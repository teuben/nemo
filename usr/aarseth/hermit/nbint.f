      SUBROUTINE NBINT(I)
*
*
*       N-body integration.
*       -------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),DV(3)
*
*
*       Obtain total force & first derivative.
      DO 5 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    5 CONTINUE
*
*       Sum over all particles.
      DO 10 J = 1,NMASS
          IF (J.EQ.I) GO TO 10
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
   10 CONTINUE
*
*       Include the corrector and set new T0, F, FDOT, D2 & D3.
      DT = TIME - T0(I)
      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
      T0(I) = TIME
*
      DO 20 K = 1,3
	  DF = 2.0*F(K,I) - FIRR(K)
	  FID = 6.0*FDOT(K,I)
	  SUM = FID + FD(K)
	  AT3 = 2.0*DF + DT*SUM
	  BT2 = -3.0*DF - DT*(SUM + FID)
*
          X0(K,I) = XI(K) + (0.6*AT3 + BT2)*DTSQ12
          X0DOT(K,I) = XIDOT(K) + (0.75*AT3 + BT2)*DT13
*
	  F(K,I) = 0.5*FIRR(K)
	  FDOT(K,I) = ONE6*FD(K)
*
          D3(K,I) = AT3*DT6
          D2(K,I) = (3.0*AT3 + BT2)*DT2
*       NOTE: These are real derivatives!
   20 continue
*
*       Specify new time-step (standard criterion or fast expression).
      IF (KZ(5).EQ.0) THEN
          TTMP = TSTEP(FIRR,FD,D2(1,I),D3(1,I),ETA)
      ELSE
          TTMP = STEPI(FIRR,FD,D2(1,I),D3(1,I),ETA)
      END IF
      DT0 = TTMP
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN 
              TTMP = MIN(2.0*STEP(I),1.0D0)
          ELSE
              TTMP = STEP(I) 
          END IF
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
          IF (TTMP.GT.DT0) THEN
              TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEP(I)
      END IF
*
*       Set new block step and update next time.
      STEP(I) = TTMP
      TNEXT(I) = STEP(I) + T0(I)
*
      RETURN
*
      END
