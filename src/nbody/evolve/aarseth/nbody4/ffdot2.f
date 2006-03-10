      SUBROUTINE FFDOT2(I,J)
*
*
*       Force & first derivative corrections.
*       -------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3)
*
*
*       Obtain dominant force and first derivative due to body #J.
      DO 15 K = 1,3
          A(K) = X(K,J) - X(K,I)
          A(K+3) = XDOT(K,J) - XDOT(K,I)
   15 CONTINUE
*
      A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
      A(8) = BODY(J)*A(7)*SQRT(A(7))
      A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
      DO 20 K = 1,3
          F1(K) = A(K)*A(8)
          F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20 CONTINUE
*
*       Subtract dominant contributions from F/2, FDOT/6 and PHI.
      DO 25 K = 1,3
          F(K,I) = F(K,I) - 0.5D0*F1(K)
          FDOT(K,I) = FDOT(K,I) - ONE6*F1DOT(K)
   25 CONTINUE
*
      PHI(I) = PHI(I) + BODY(J)*SQRT(A(7))
*
      RETURN
*
      END
