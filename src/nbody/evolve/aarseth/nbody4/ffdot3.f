      SUBROUTINE FFDOT3(I,KCASE)
*
*
*       Force & first derivative from neighbours.
*       -----------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3)
*
*
*       Initialize force & first derivative (KCASE = 2) of body #I.
      DO 10 K = 1,3
          IF (KCASE.EQ.2) F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all neighbours.
      NNB = ILIST(1)
      DO 30 L = 2,NNB+1
          J = ILIST(L)
          IF (J.EQ.I) GO TO 30
*
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2
          A(7) = 1.0/RIJ2
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          DO 20 K = 1,3
              IF (KCASE.EQ.2) F1(K) = A(K)*A(8)
              F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
*
          DO 25 K = 1,3
              IF (KCASE.EQ.2) F(K,I) = F(K,I) + F1(K)
              FDOT(K,I) = FDOT(K,I) + F1DOT(K)
   25     CONTINUE
*
*       Include close neighbours inside 5*RMIN in the potential.
          IF (KCASE.EQ.2.AND.RIJ2.LT.25.0*RMIN2) THEN
              PHI(I) = PHI(I) - BODY(J)/SQRT(RIJ2)
          END IF
*
   30 CONTINUE
*
      RETURN
*
      END
