      SUBROUTINE FPOLY1(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common1.h'
      REAL*4  A(6)
*
*
*       Loop over all bodies or just one single body (I1 = I2).
      DO 50 I = I1,I2
*
*       Initialize force & first derivative.
      DO 2 K = 1,3
          F(K,I) = 0.0
          FDOT(K,I) = 0.0
    2 CONTINUE
*
*       Sum over all the other particles.
      DO 10 J = 1,N
          IF (J.EQ.I) GO TO 10
*
          DO 5 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
    5     CONTINUE
*
          A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3) + EPS2)
          A8 = BODY(J)*A7*SQRT(A7)
          A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A7
*
          DO 8 K = 1,3
              F(K,I) = F(K,I) + A(K)*A8
              FDOT(K,I) = FDOT(K,I) + (A(K+3) - A(K)*A9)*A8
    8     CONTINUE
   10 CONTINUE
   50 CONTINUE
*
      RETURN
*
      END
