      SUBROUTINE NBINT(I,IR,XI,YI,ZI,XIDOT)
*
*
*       Irregular integration.
*       ----------------------
*
      INCLUDE 'common2.h'
      REAL*8  XI,YI,ZI,FIRR(3)
      REAL*4  XIDOT(3)
*
*
*       Set current neighbour number.
      NNB0 = LIST(1,I)
      NNB1 = NNB0 + 1
*
*       See whether the regular force needs to be updated.
      IF (TIME + STEP(I).GT.T0R(I) + STEPR(I)) GO TO 40 
*
*       Predict the coordinates of neighbours to order FDOT.
      IR = 0
   10 DO 20 L = 2,NNB1
          J = LIST(L,I)
          S = TIME - T0(J)
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + XDOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + XDOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + XDOT(3,J))*S + X0(3,J)
   20 CONTINUE
*
      GO TO 100
*
*       Set non-zero indicator for updating regular force field.
   40 IR = 1
*       See whether to predict all particles or just the neighbours.
      IF (KZ(14).GT.0.AND.NNB0.GT.KZ(14)) GO TO 10
*
*       Only perform full coordinate prediction for small NNB if KZ(14) > 0.
      NNPRED = NNPRED + 1
*
      DO 50 J = 1,N
          S = TIME - T0(J)
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + XDOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + XDOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + XDOT(3,J))*S + X0(3,J)
   50 CONTINUE
*
*       Predict coordinates and velocities of body #I to order F3DOT.
  100 DT = TIME - T0(I)
      T1PR = T0(I) - T1(I)
      T2PR = T0(I) - T2(I)
      T12PR = T1PR + T2PR
      DT06 = 0.6*DT
      DT19 = ONE9*DT
      DT12 = ONE12*DT
      DT34 = 0.75*DT
      DT32 = 1.5*DT
      DT20 = 2.0*DT  
      DT0R = (T0(I) - T0R(I)) + (T0(I) - T1R(I)) + (T0(I) - T2R(I))
*
      DO 115 K = 1,3
          F2DOTK = D3R(K,I)*DT0R + D2R(K,I) + D3(K,I)*T12PR + D2(K,I)
          F3DOTK = D3R(K,I) + D3(K,I)
          X(K,I) = ((((F3DOTK*DT06 + F2DOTK)*DT12 + FDOT(K,I))*DT + 
     &                             F(K,I))*DT + X0DOT(K,I))*DT + X0(K,I)
          X0DOT(K,I) = (((F3DOTK*DT34 + F2DOTK)*DT19 + FDOT(K,I))*DT32 +
     &                                         F(K,I))*DT20 + X0DOT(K,I) 
*       Save velocity without semi-iteration for derivative corrections.
          XIDOT(K) = X0DOT(K,I)
          FIRR(K) = 0.0
  115 CONTINUE
*
*       Obtain irregular force.
      XI = X(1,I)
      YI = X(2,I)
      ZI = X(3,I)
*
*       Assume small mass at centre for special case of no neighbours.
      IF (NNB0.EQ.0) THEN
          RI2 = (XI - RDENS(1))**2 + (YI - RDENS(2))**2 +
     &                               (ZI - RDENS(3))**2
          A6 = 0.01*BODYM/(RI2*SQRT(RI2))
          DO 122 K = 1,3
              FIRR(K) = -A6*X(K,I)
  122     CONTINUE
          GO TO 170
      END IF
*
*       Set maximum close encounter distance.
      RCL2 = 5.01*EPS2 
*
*       Sum all neighbour contributions.
      DO 140 L = 2,NNB1
          K = LIST(L,I)
          A1 = X(1,K) - XI
          A2 = X(2,K) - YI
          A3 = X(3,K) - ZI
          RIJ2 = A1*A1 + A2*A2 + A3*A3 + EPS2
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
*
*       ------------------------------------------------
*       Check possible merger candidate (suppress if KZ(12) not used).
          IF (RIJ2.LT.RCL2) THEN
              RCL2 = RIJ2
              JCOMP = K
          END IF
*       ------------------------------------------------
*
*       Accumulate the force in double precision to improve accuracy.
          FIRR(1) = FIRR(1) + A1*A6
          FIRR(2) = FIRR(2) + A2*A6
          FIRR(3) = FIRR(3) + A3*A6
  140 CONTINUE
*
*       ------------------------------------------------
*       Check merger candidate (RIJ < 2*SQRT(EPS2).
      IF (KZ(12).GT.0.AND.RCL2 - EPS2.LT.4.0*EPS2) THEN
          NMTRY = NMTRY + 1
*       Obtain semi-major axis.
          VIJ2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 +
     &           (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &           (XDOT(3,I) - XDOT(3,JCOMP))**2
          SEMI = 2.0/SQRT(RCL2) - VIJ2/(BODY(I) + BODY(JCOMP))
          SEMI = 1.0/SEMI
*       Adopt inelastic merger for significantly bound two-body motion.
          IF (SEMI.GT.0.0.AND.SEMI.LT.0.05*RSCALE) THEN
*       Save primary body and set merger indicator for end of cycle.
              ICOMP = I
              IPHASE = 4
          END IF
      END IF
*       ------------------------------------------------
*
*       Set time intervals for corrector and update irregular force times.
  170 DT1 = TIME - T1(I)
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
*       Form new differences and include semi-iteration for irregular force.
      DO 180 K = 1,3
          D1K = (FIRR(K) - FI(K,I))*A1
          D2K = (D1K - D1(K,I))*A2
          D3K = (D2K - D2(K,I))*A3
          F4DOTK = (D3K - D3(K,I))*A4
          FI(K,I) = FIRR(K)
          D1(K,I) = D1K
          D2(K,I) = D2K
          D3(K,I) = D3K
          X(K,I) = F4DOTK*S6 + X(K,I)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = F4DOTK*S7 + X0DOT(K,I)
          XDOT(K,I) = X0DOT(K,I)
  180 CONTINUE
*
*       See whether a new regular force is required.
      IF (IR.EQ.0) THEN
*       Form time-step factors for prediction variables.
          DTR = TIME - T0R(I)
          DT1R = TIME - T1R(I)
          DT2R = TIME - T2R(I)
          S1 = DTR + DT1R
          S2 = DT2R*S1 + DTR*DT1R
*
*       Extrapolate regular force & first derivatives to obtain F & FDOT.
          DO 190 K = 1,3
              F1K = ((D3R(K,I)*DT2R + D2R(K,I))*DT1R + D1R(K,I))*DTR +
     &                                                           FR(K,I)
              F1DOTK = D3R(K,I)*S2 + D2R(K,I)*S1 + D1R(K,I) +
     &                              (D3(K,I)*DT1 + D2(K,I))*DT + D1(K,I)
              F(K,I) = 0.5*(F1K + FIRR(K))
              FDOT(K,I) = ONE6*F1DOTK
  190     CONTINUE
      END IF
*
      NSTEPI = NSTEPI + 1
*
      RETURN
*
      END
