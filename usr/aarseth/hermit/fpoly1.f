      SUBROUTINE FPOLY1(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  A(9),F1(3),F1DOT(3),XI(3),XIDOT(3)
*
*
*       Loop over all bodies or one single body.
      DO 40 I = I1,I2
*
*       Initialize force & first derivative for body #I.
      DO 10 K = 1,3
          F1(K) = 0.0
          F1DOT(K) = 0.0
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all bodies.
      DO 30 J = 1,N
          IF (J.EQ.I) GO TO 30
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          DO 20 K = 1,3
              F1(K) = F1(K) + A(K)*A(8)
              F1DOT(K) = F1DOT(K) + (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
   30 CONTINUE
*
      DO 35 K = 1,3
          F(K,I) = F1(K)
          FDOT(K,I) = F1DOT(K)
   35 CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
