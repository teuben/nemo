      SUBROUTINE XVPRED(I1,NNB)
*
*
*       Prediction of coordinates & velocities.
*       ---------------------------------------
*
      INCLUDE 'commonp.h'
*
*
*       Perform full N prediction to highest order (defined by NNB >= N).
      IF (NNB.GE.N) GO TO 20
*       Use high order for special singles (NNB < 0).
      IF (NNB.LT.0) GO TO 20
*
*       Adopt low order for standard single particles (I1 <= N & NNB = 0).
      IF (NNB.EQ.0) THEN
          I = I1
          S = TIME - T0(I)
          IF (S.EQ.0.0D0) GO TO 40
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
          X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
          X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
          XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
          XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
          XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
          GO TO 40
      END IF
*
*       Set index of first body (#I1 if NNB = 0 or NNB >= N).
   20 I = I1
*       Predict coordinates & velocities of body #I to order F3DOT.
   25 DT = TIME - T0(I)
      IF (DT.EQ.0.0D0) GO TO 35
      A1 = 0.5*DT
      A2 = ONE12*DT
      A3 = ONE6*DT
      A4 = 0.25*DT
*
      DO 30 K = 1,3
          X(K,I) = ((((D3(K,I)*A1 + D2(K,I))*A2 + FDOT(K,I))*DT +
     &                 F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
          XDOT(K,I)  = (((D3(K,I)*A4 + D2(K,I))*A3 + 3.0*FDOT(K,I))*DT
     &                               + 2.0*F(K,I))*DT + X0DOT(K,I)
   30 CONTINUE
*
*       Check whether prediction of single particle or full N.
   35 IF (NNB.GT.0) THEN
          I = I + 1
          IF (I.LE.NNB) GO TO 25
      END IF
*
   40 RETURN
*
      END
