      SUBROUTINE NBINT(I)
*
*
*       N-body integrator.
*       ------------------
*
      INCLUDE 'common1.h'
      REAL*8  XI,YI,ZI
      REAL*4  FIRR(3),F1DOT(4)
*
*
*       Predict all coordinates to first order in force derivative.
      DO 7 J = 1,N
          S = TIME - T0(J)
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + XDOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + XDOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + XDOT(3,J))*S + X0(3,J)
    7 CONTINUE
*
*       Improve coordinates & velocities of body #I to order F3DOT.
      DT = TIME - T0(I)
      T1PR = T0(I) - T1(I)
      T2PR = T0(I) - T2(I)
      T12PR = T1PR + T2PR
      DT06 = 0.6*DT
      DT19 = ONE9*DT
      DT12 = ONE12*DT
      DT34 = 0.75*DT
      DT32 = 1.5*DT
      DT20 = 2.0*DT  
*
      DO 9 K = 1,3
          F2DOTK = D3(K,I)*T12PR + D2(K,I)
          X(K,I) = ((((D3(K,I)*DT06 + F2DOTK)*DT12 + FDOT(K,I))*DT + 
     &                             F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
          X0DOT(K,I) = (((D3(K,I)*DT34 + F2DOTK)*DT19 +
     &                       FDOT(K,I))*DT32 + F(K,I))*DT20 + X0DOT(K,I) 
          FIRR(K) = 0.0
    9 CONTINUE
*
*       Obtain the current force.
      XI = X(1,I)
      YI = X(2,I)
      ZI = X(3,I)
*
      DO 10 J = 1,N
          IF (J.EQ.I) GO TO 10
          A1 = X(1,J) - XI
          A2 = X(2,J) - YI
          A3 = X(3,J) - ZI
          RIJ2 = A1*A1 + A2*A2 + A3*A3 + EPS2
          A5 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FIRR(1) = FIRR(1) + A1*A5
          FIRR(2) = FIRR(2) + A2*A5
          FIRR(3) = FIRR(3) + A3*A5
   10 CONTINUE
*
*       Set time intervals for corrector and update force times.
      DT1 = TIME - T1(I)
      DT2 = TIME - T2(I)
      DT3 = TIME - T3(I)
      T3PR = T0(I) - T3(I)
      S2 = T1PR*T2PR
      S3 = S2*T3PR
      S4 = S2 + T3PR*T12PR
      S5 = T12PR + T3PR
      S6 = (((0.6666667*DT + S5)*DT06 + S4)*DT12 + ONE6*S3)*DT
      S7 = ((0.2*DT + 0.25*S5)*DT + ONE3*S4)*DT + 0.5*S3
      T3(I) = T2(I)
      T2(I) = T1(I)
      T1(I) = T0(I)
      T0(I) = TIME
      A1 = 1.0/DT
      A2 = 1.0/DT1
      A3 = 1.0/DT2
      A4 = DT*DT/DT3
*
*       Form new differences and include fourth-order semi-iteration.
      DO 20 K = 1,3
          D1K = (FIRR(K) - 2.0*F(K,I))*A1
          D2K = (D1K - D1(K,I))*A2
          D3K = (D2K - D2(K,I))*A3
          F4DOTK = (D3K - D3(K,I))*A4
          D1(K,I) = D1K 
          D2(K,I) = D2K 
          D3(K,I) = D3K 
          X(K,I) = F4DOTK*S6 + X(K,I)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = F4DOTK*S7 + X0DOT(K,I)
          XDOT(K,I) = X0DOT(K,I)
   20 CONTINUE
*
*       Set half the force & sixth the first derivative for fast prediction.
      DO 30 K = 1,3
          F1DOT(K) = (D3(K,I)*DT1 + D2(K,I))*DT + D1(K,I)
          F(K,I) = 0.5*FIRR(K)
          FDOT(K,I) = ONE6*F1DOT(K)
   30 CONTINUE
*
      NSTEPI = NSTEPI + 1
*
      RETURN
*
      END
