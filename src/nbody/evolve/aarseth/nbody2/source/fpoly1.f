      SUBROUTINE FPOLY1(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(6),F1DOT(3)
      REAL*8  F1(3)
*
*
*       Loop over all bodies or just one single body (I1 = I2).
      DO 40 I = I1,I2
*
*       Initialize irregular & regular force and first derivative.
      DO 10 K = 1,3
          FI(K,I) = 0.0
          D1(K,I) = 0.0
          FR(K,I) = 0.0
          D1R(K,I) = 0.0
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all bodies.
      NNB = LIST(1,I)
      L = 2
      NAMEJ = LIST(L,I)
*       Index of first neighbour to be identified in force loop.
*
      DO 30 J = 1,N
          IF (J.EQ.I) GO TO 30
*
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A7 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3) + EPS2) 
          A8 = BODY(J)*A7*SQRT(A7)
          A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A7
*
          DO 20 K = 1,3
              F1(K) = A(K)*A8
              F1DOT(K) = (A(K+3) - A(K)*A9)*A8
   20     CONTINUE
*
          IF (J.NE.NAMEJ) THEN
              DO 25 K = 1,3
                  FR(K,I) = FR(K,I) + F1(K)
                  D1R(K,I) = D1R(K,I) + F1DOT(K)
   25         CONTINUE
          ELSE
              DO 28 K = 1,3
                  FI(K,I) = FI(K,I) + F1(K)
                  D1(K,I) = D1(K,I) + F1DOT(K)
   28         CONTINUE
*       Advance the neighbour list until last member has been considered.
              IF (L.LE.NNB) THEN
                  L = L + 1
                  NAMEJ = LIST(L,I)
              END IF
          END IF
   30 CONTINUE
*
*       Check external force option (F1 is dummy argument).
      IF (KZ(15).GT.0) THEN
          IF (KZ(15).EQ.1) THEN
              CALL XTRNL1(I,2,F1)
          ELSE
              CALL XTRNL2(I,2,F1)
          END IF 
      END IF
*
*       Set total force & first derivative.
      DO 35 K = 1,3
          F(K,I) = FI(K,I) + FR(K,I)
          FDOT(K,I) = D1(K,I) + D1R(K,I)
   35 CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
