      SUBROUTINE INTGRT
*
*
      INCLUDE 'common1.h'
      REAL*8  k0,k1,k2,k3,k5
      SAVE
*
*
*       Skip the first part at later times.
      IF (TIME.EQ.0.0D0) THEN
*       Initialize previous step and specify small initial step.
          h1 = STEP
          h2 = 1.2*STEP
          STEP0 = h1
      END IF
*
*       Form integration coefficients.
  10  k0 = 0.5*(h1 + h2)
      k1 = h2**3/24.0D0
      k2 = (h2**2 - h1**2)/8.0D0
      k3 = (h2**2 + 2.0*h1*h2 - 2.0*h1**2)/24.0D0
*       Suppress the original bad expression (cf. Z.f.Astrophys).
*     k3 = (h2**2 + 2.0*h2*h1 - h1**2)/24.0D0
      k5 = h2*(h2**2 + h2*h1 - h1**2)/12.0D0
*
*       Predict coordinates and velocities to order FDOT.
      DO 25 I = 1,N
          DO 20 K = 1,3
              F1(K,I) = b1(K,I)
              b1(K,I) = F(K,I)
              a1(K,I) = (b1(K,I) - F1(K,I))/h1
              XDOT(K,I) = XDOT(K,I) + k0*b1(K,I) + k2*a1(K,I)
              X(K,I) = X(K,I) + h2*XDOT(K,I) + k1*a1(K,I)
   20     CONTINUE
   25 CONTINUE
*
*       Obtain forces at end of the step.
      CALL FORCE
*
*       Evaluate new force differences and apply the corrector.
      DO 40 I = 1,N
          DO 30 K = 1,3
              b2(K,I) = F(K,I)
              a2(K,I) = (b2(K,I) - b1(K,I))/h2
              d2(K,I) = a2(K,I) - a1(K,I)
              XDOT(K,I) = XDOT(K,I) + k3*d2(K,I)
              X(K,I) = X(K,I) + k5*d2(K,I)
   30     CONTINUE
   40 CONTINUE
*
*       Advance the time and update time-steps.
      TIME = TIME + h2
      STEP1 = STEP0
      STEP0 = h2
      h1 = h2
      h2 = MIN(STEP,1.2*h2)
      NSTEPI = NSTEPI + 1
*
      IF (TIME.LT.TNEXT) GO TO 10
*
      RETURN
*
      END
