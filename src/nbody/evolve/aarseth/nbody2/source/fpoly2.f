      SUBROUTINE FPOLY2(I1,I2)
*
*
*       Second & third force derivative.
*       --------------------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(12),F2DOT(3),F3DOT(3)
*
*
*       Loop over all bodies or just one single body (I1 = I2).
      DO 70 I = I1,I2
*
*       Obtain second & third force derivative. 
      NNB = LIST(1,I)
*
*       Neglect F2DOT & F3DOT outside 3*RS unless high accuracy is needed.
      RCRIT2 = 9.0*RS(I)**2*(1.0 + 1.0/FLOAT(NNB))
*
*       Include all contributions if there is no unique density centre.
      IF (KZ(10).GT.0) THEN
          RCRIT2 = (100.0*RS(I))**2
      END IF
*
*       Initialize the higher differences for body #I.
      DO 10 K = 1,3
          D2(K,I) = 0.0
          D3(K,I) = 0.0
          D2R(K,I) = 0.0
          D3R(K,I) = 0.0
   10 CONTINUE
      L = 2
      NAMEJ = LIST(L,I)
*
      DO 60 J = 1,N 
          IF (J.EQ.I) GO TO 60
*
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
   15     CONTINUE
          RIJ2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
*
*       Ensure that all neighbours are considered.
          IF (RIJ2.GT.RCRIT2.AND.J.NE.NAMEJ) GO TO 60
*       Distant bodies do not contribute significantly to F2DOT & F3DOT.
*
          IF (TIME.LE.0.0D0) THEN
*       Copy current force & derivative if TIME = 0.
              DO 20 K = 1,3
                  A(K+6) = F(K,J)
                  A(K+9) = FDOT(K,J)
   20         CONTINUE
          ELSE
*       Obtain current force and first derivative to second order.
              DT = TIME - T0(J)
              DT1 = TIME - T1(J)
              DTR = TIME - T0R(J)
              DT1R = TIME - T1R(J)
*
              DO 35 K = 1,3
                  A(K+6) = (D2R(K,J)*DT1R + D1R(K,J))*DTR + FR(K,J) +
     &                              (D2(K,J)*DT1 + D1(K,J))*DT + FI(K,J)
                  A(K+9) = D2R(K,J)*(DTR + DT1R) + D1R(K,J) +
     &                                      D2(K,J)*(DT + DT1) + D1(K,J)
   35         CONTINUE
          END IF
*
          DO 45 K = 1,3
              A(K+3) = XDOT(K,J) - XDOT(K,I)
              A(K+6) = A(K+6) - F(K,I)
              A(K+9) = A(K+9) - FDOT(K,I)
   45     CONTINUE
*
          A13 = 1.0/(RIJ2 + EPS2)
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
   50     CONTINUE
*
*       See whether body #J is a neighbour of body #I.
          IF (J.NE.NAMEJ) THEN
              DO 52 K = 1,3
                  D2R(K,I) = D2R(K,I) + F2DOT(K)
                  D3R(K,I) = D3R(K,I) + F3DOT(K)
   52         CONTINUE
          ELSE
              DO 55 K = 1,3
                  D2(K,I) = D2(K,I) + F2DOT(K)
                  D3(K,I) = D3(K,I) + F3DOT(K)
   55         CONTINUE
*
*       Advance the neighbour list until last member is identified.
              IF (L.LE.NNB) THEN
                  L = L + 1
                  NAMEJ = LIST(L,I)
              END IF
          END IF
   60 CONTINUE
   70 CONTINUE
*
      RETURN
*
      END
