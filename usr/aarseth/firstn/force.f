      SUBROUTINE FORCE
*
*
*       Total force and time-step.
*       --------------------------
*
      INCLUDE 'common1.h'
      REAL*4  A(3)
*
*
*       Loop over all bodies.
      RM2 = 100.0
      DO 50 I = 1,N
*
*       Initialize the force.
      DO 2 K = 1,3
          F(K,I) = 0.0
    2 CONTINUE
*
*       Sum over all the other particles.
      DO 10 J = 1,N
          IF (J.EQ.I) GO TO 10
*
          DO 5 K = 1,3
              A(K) = X(K,J) - X(K,I)
    5     CONTINUE
*
          RIJ2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3) + EPS2
*       Find the closest particle pair.
          IF (RIJ2.LT.RM2) THEN
              RM2 = RIJ2
              ICL = I
              JCL = J
          END IF
          A8 = BODY(J)/(RIJ2*SQRT(RIJ2))
*
          DO 8 K = 1,3
              F(K,I) = F(K,I) + A(K)*A8
    8     CONTINUE
   10 CONTINUE
   50 CONTINUE
*
*       Set the time-step from smallest separation and relative velocity.
      VIJ2 = 0.0
      DO 60 K = 1,3
          VIJ2 = VIJ2 + (XDOT(K,ICL) - XDOT(K,JCL))**2
   60 CONTINUE
      RM = SQRT(RM2 - EPS2)
      VM = SQRT(VIJ2)
      STEP = RM*SQRT(2.0*RM)/(1.0 + VM*SQRT(2.0*RM))
      STEP = STEP/ETA
*
*       Specify smaller initial step.
*     IF (TIME.EQ.0.0D0) THEN
*         STEP = 0.5*STEP
*         STEP1 = 0.8*STEP
*         STEP0 = STEP
*     END IF
*
      RETURN
*
      END
