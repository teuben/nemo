      SUBROUTINE STEPS(I1,I2)
*
*
*       Initialization of time-steps & differences.
*       -------------------------------------------
*
      INCLUDE 'common2.h'
      REAL*4  F2DOT(4),F3DOT(4)
*
*
*       Set new steps and initialize divided differences.
      DO 40 I = I1,I2
*
*       Obtain irregular & regular time-steps using the composite formula.
      FI2 = FI(1,I)**2 + FI(2,I)**2 + FI(3,I)**2
      DI2 = D1(1,I)**2 + D1(2,I)**2 + D1(3,I)**2
*
      DO 10 K = 1,3
          F2DOT(K) = D2(K,I)
          F3DOT(K) = D3(K,I)
   10 CONTINUE
*
      F2DOT(4) = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      F3DOT(4) = F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2
      A1 = (SQRT(FI2*F2DOT(4)) + DI2)/
     &                             (SQRT(DI2)*SQRT(F3DOT(4)) + F2DOT(4))
      STEP(I) = SQRT(ETAI*A1)
*
      FR2 = FR(1,I)**2 + FR(2,I)**2 + FR(3,I)**2
      DR2 = D1R(1,I)**2 + D1R(2,I)**2 + D1R(3,I)**2
*
      DO 20 K = 1,3
          F2DOT(K) = D2R(K,I)
          F3DOT(K) = D3R(K,I)
   20 CONTINUE
*
      F2DOT(4) = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      F3DOT(4) = F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2
*
*       Use dominant second derivative (F2DOT = FR**2/RS) in rare case.
      IF (F2DOT(4).LT.(FR2/RS(I))**2) THEN
          F2DOT(4) = (FR2/RS(I))**2
      END IF 
*
      A1 = (SQRT(FR2*F2DOT(4)) + DR2)/
     &                             (SQRT(DR2)*SQRT(F3DOT(4)) + F2DOT(4))
      STEPR(I) = SQRT(ETAR*A1)
*
*       Ensure that regular step is not too large (DTR < 0.1*RS/V).
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      IF (VI2.GT.0.0) THEN
          STEPR(I) = MIN(STEPR(I),0.1*RS(I)/SQRT(VI2))
      END IF
*       Perform safety test on irregular step.
      STEP(I) = MIN(STEPR(I),STEP(I))
*
*       Initialize backward times.
      T0(I) = TIME
      T1(I) = TIME - STEP(I)
      T2(I) = TIME - 2.0*STEP(I)
      T3(I) = TIME - 3.0*STEP(I)
      T0R(I) = TIME
      T1R(I) = TIME - STEPR(I)
      T2R(I) = TIME - 2.0*STEPR(I)
      T3R(I) = TIME - 3.0*STEPR(I)
      DT1 = STEP(I)
      DT1R = STEPR(I)
*
*       Convert from derivatives to divided differences.
      DO 30 K = 1,3
          D1(K,I) = (ONE6*D3(K,I)*DT1 - 0.5*D2(K,I))*DT1 + D1(K,I)
          D2(K,I) = 0.5*D2(K,I) - 0.5*D3(K,I)*DT1
          D3(K,I) = ONE6*D3(K,I)
          D1R(K,I) = (ONE6*D3R(K,I)*DT1R - 0.5*D2R(K,I))*DT1R + D1R(K,I)
          D2R(K,I) = 0.5*D2R(K,I) - 0.5*D3R(K,I)*DT1R
          D3R(K,I) = ONE6*D3R(K,I)
*
*       Initialize the primary coordinates & velocities.
          X0(K,I) = X(K,I)
          X0DOT(K,I) = XDOT(K,I)
*
*       Set half the force & sixth the first derivative for fast prediction.
          F(K,I) = 0.5*F(K,I)
          FDOT(K,I) = ONE6*FDOT(K,I)
   30 CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
