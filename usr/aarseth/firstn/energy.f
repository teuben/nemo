      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common1.h'
      COMMON/VELOC/  X2DOT(3,NMAX)
      REAL*8  XD(3)
*
*
*       Calculate the potential energy.
      ZKIN = 0.0
      POT = 0.0
      I = 1
   10 POTJ = 0.0
      DO 15 J = I+1,N
          X1 = X(1,I) - X(1,J)
          X2 = X(2,I) - X(2,J)
          X3 = X(3,I) - X(3,J)
          POTJ = POTJ + BODY(J)/SQRT(X1*X1 + X2*X2 + X3*X3 + EPS2)
   15 CONTINUE
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.N) GO TO 10
*
*       Assign appropriate steps (NB! two previous values).
      IF (TIME.EQ.0.0D0) THEN
          STEP1 = STEP
          STEP0 = STEP
      END IF
      h1 = STEP1
      h2 = STEP0
*
*       Obtain the kinetic energy (using average force is not sufficient).
      DO 30 I = 1,N
          DO 25 K = 1,3
              IF (TIME.EQ.0.0D0) b1(K,I) = F(K,I)
*             XD(K) = XDOT(K,I) + 0.25*(b1(K,I) + F(K,I))*DT
              XD(K) = XDOT(K,I) + 0.5*h2*b1(K,I) + 3.0*h2**2*a1(K,I)/8.0
     &                + d2(K,I)*h2**2*(3.0*h1/8.0 + 7.0*h2/24.0)/(h1+h2)
              ZKIN = ZKIN + BODY(I)*XD(K)**2
*       Copy current velocity for accuracy checks.
              X2DOT(K,I) = XD(K)
   25     CONTINUE
   30 CONTINUE
*
      ZKIN = 0.5*ZKIN
*       Total energy = ZKIN - POT.
*
      RETURN
*
      END
