      SUBROUTINE FPOLY2(I1,I2)
*
*
*       Second & third force derivative.
*       --------------------------------
*
      INCLUDE 'common1.h'
      REAL*4  A(12),F2DOT(3),F3DOT(3)
*
*
*       Loop over all bodies or just one single body (I1 = I2).
      DO 70 I = I1,I2
*
*       Initialize second & third derivative.
      DO 10 K = 1,3
          D2(K,I) = 0.0
          D3(K,I) = 0.0
   10 CONTINUE
*
*       Obtain second and third force derivative.
      DO 60 J = 1,N
          IF (J.EQ.I) GO TO 60
*
          DO 35 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
              A(K+6) = F(K,J) - F(K,I)
              A(K+9) = FDOT(K,J) - FDOT(K,I)
   35     CONTINUE
*
          A13 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3) + EPS2)
          A14 = BODY(J)*A13*SQRT(A13)
          A15 = (A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A13
          A16 = A15*A15
          A17 = 3.0*A15
          A18 = 6.0*A15
          A19 = 9.0*A15
          A20 = (A(4)*A(4) + A(5)*A(5) + A(6)*A(6) + A(1)*A(7) +
     &                                  A(2)*A(8) + A(3)*A(9))*A13 + A16
          A21 = 9.0*A20
          A20 = 3.0*A20
          A22 = (9.0*(A(4)*A(7) + A(5)*A(8) + A(6)*A(9)) +
     &                 3.0*(A(1)*A(10) + A(2)*A(11) + A(3)*A(12)))*A13 +
     &                                               A17*(A20 - 4.0*A16)
*
          DO 50 K = 1,3
              F1DOTK = A(K+3) - A17*A(K)
              F2DOT(K) = (A(K+6) - A18*F1DOTK - A20*A(K))*A14
              F3DOT(K) = (A(K+9) - A21*F1DOTK - A22*A(K))*A14 -
     &                                                      A19*F2DOT(K)
              D2(K,I) = D2(K,I) + F2DOT(K)
              D3(K,I) = D3(K,I) + F3DOT(K)
   50     CONTINUE
   60 CONTINUE
   70 CONTINUE
*
      RETURN
*
      END
