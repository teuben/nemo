      SUBROUTINE FFDOT(I)
*
*
*       Force & first derivative from neighbours.
*       -----------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3)
*
*
*       Initialize force & first derivative of body #I.
      DO 10 K = 1,3
          F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all perturbers.
      KDUM = 0
      NNB = ILIST(1)
      DO 30 L = 2,NNB+1
          J = ILIST(L)
          IF (J.EQ.I) GO TO 30
          IF (J.GT.N) THEN
              JPAIR = J - N
*       Resolve perturbed KS or use c.m. approximation for unperturbed case.
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  ZZ = 1.0
                  IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0D0
                  CALL KSRES2(JPAIR,J1,J2,ZZ)
                  KDUM = J1
                  J = KDUM
              END IF
          END IF
*
   12     DO 15 K = 1,3
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
              F1(K) = A(K)*A(8)
              F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
*
          DO 25 K = 1,3
              F(K,I) = F(K,I) + F1(K)
              FDOT(K,I) = FDOT(K,I) + F1DOT(K)
   25     CONTINUE
*
*       Include close encounter effect in the potential.
          IF (RIJ2.LT.25.0*RMIN2) THEN
              PHI(I) = PHI(I) - BODY(J)/SQRT(RIJ2)
          END IF
*
          IF (J.EQ.KDUM) THEN
              J = J + 1
              GO TO 12
          END IF
   30 CONTINUE
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(I,I,1)
      END IF
*
      RETURN
*
      END
