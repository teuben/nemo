      SUBROUTINE STEPS(I1,I2)
*
*
*       Initialization of time-steps & differences.
*       -------------------------------------------
*
      INCLUDE 'common1.h'
      REAL*4  F2DOT(4),F3DOT(4)
*
*
*       Set new steps and initialize divided differences.
      DO 40 I = I1,I2
*
*       Initialize integration steps and convert to force differences.
      FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
      DI2 = FDOT(1,I)**2 + FDOT(2,I)**2 + FDOT(3,I)**2
*
      DO 10 K = 1,3
          F2DOT(K) = D2(K,I)
          F3DOT(K) = D3(K,I)
   10 CONTINUE
*
      F2DOT(4) = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      F3DOT(4) = F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2
      DT2 = (SQRT(FI2*F2DOT(4)) + DI2)/(SQRT(DI2*F3DOT(4)) + F2DOT(4))
      STEP(I) = SQRT(ETA*DT2)
*
*       Initialize the backward times.
      T0(I) = TIME
      T1(I) = TIME - STEP(I)
      T2(I) = TIME - 2.0*STEP(I)
      T3(I) = TIME - 3.0*STEP(I)
      DT1 = STEP(I)
*
*       Convert from derivatives to divided differences.
      DO 30 K = 1,3
          D1(K,I) = (ONE6*D3(K,I)*DT1 - 0.5*D2(K,I))*DT1 + FDOT(K,I)
          D2(K,I) = 0.5*D2(K,I) - 0.5*D3(K,I)*DT1
          D3(K,I) = ONE6*D3(K,I)
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
