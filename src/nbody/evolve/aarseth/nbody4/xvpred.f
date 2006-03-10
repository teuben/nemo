      SUBROUTINE XVPRED(I1,I2)
*
*
*       Prediction of coordinates & velocities.
*       ---------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Choose between low or high-order prediction (I2 = -2; I2 >= N).
      IF (I2.EQ.-2) GO TO 1
      IF (I2.GE.N) GO TO 20
*       Exclude special cases (I2 < 0: single particle; I1 > N: c.m. body).
      IF (I2.LT.0.OR.I1.GT.N) GO TO 20
*
*       Predict sequential bodies or just one particle (I2 >= I1; I2 = 0).
    1 I3 = MAX(I1,I2)
      DO 5 I = I1,I3
          S = TIME - T0(I)
          IF (S.EQ.0.0D0) GO TO 5
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
          X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
          X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
          XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
          XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
          XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
    5 CONTINUE
*
*       Resolve the components of any perturbed pairs (skip I2 < 0).
      I = I2
   10 IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
              CALL RESOLV(JPAIR,1)
          END IF
          I = I - 1
          GO TO 10
      END IF
      GO TO 40
*
*       Set index of first body.
   20 I = I1
*       Predict coordinates & velocities of body #I to order F3DOT.
*       Note: D2 now stored as D2/18 for GRAPE-6.
   25 DT = TIME - T0(I)
      IF ((DT.EQ.0.0D0.AND.IPHASE.LT.4).OR.BODY(I).EQ.0.0D0) GO TO 35
      A1 = 0.2*DT
      A2 = 0.04166666666667D0*DT
      DO 30 K = 1,3
          X(K,I) = ((((D3(K,I)*A1 + 18.D0*D2(K,I))*A2 + FDOT(K,I))*DT +
     &                             F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
          XDOT(K,I)  = (((D3(K,I)*A2 + 3.D0*D2(K,I))*DT +
     &                              3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT +
     &                                                        X0DOT(K,I)
   30 CONTINUE
*
*       Resolve the components of perturbed pairs.
   35 IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
              CALL RESOLV(JPAIR,1)
          END IF
      END IF
*
*       Check whether more particles need predicting (I2 is last if > 0).
      IF (I2.GT.0) THEN
          I = I + 1
          IF (I.LE.I2) GO TO 25
      END IF
*
   40 RETURN
*
      END
