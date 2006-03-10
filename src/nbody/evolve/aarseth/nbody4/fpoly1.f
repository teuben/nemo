      SUBROUTINE FPOLY1(I1,I2,KCASE)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3)
*
*
*       Standard case, new c.m. or KS termination (KCASE = 0, 1, 2).
      JLAST = NTOT
*       Reduce loop size for new c.m. polynomial.
      IF (KCASE.EQ.1) JLAST = NTOT - 1
*
*       Loop over all bodies, pair #ICOMP & JCOMP or one single body.
      DO 40 I = I1,I2
*
*       Initialize force & first derivative of body #I.
      DO 10 K = 1,3
          F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
   10 CONTINUE
      PHI(I) = 0.0D0
*
*       Obtain force & first derivative by summing over all bodies.
      KDUM = 0
*
      DO 30 JDUM = IFIRST,JLAST
          IF (JDUM.EQ.I) GO TO 30
          J = JDUM
          IF (J.GT.N) THEN
              JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  KDUM = 2*JPAIR - 1
                  J = KDUM
              END IF
          END IF
*
   12     DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          DO 20 K = 1,3
              F1(K) = A(K)*A(8)
              F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
*
*       See whether summation index is equal to either component.
          IF (J.EQ.ICOMP.OR.J.EQ.JCOMP) THEN
              IF (KCASE.EQ.1) GO TO 30
*       Note that dominant terms cancel analytically in c.m. polynomial.
          END IF
*
          DO 25 K = 1,3
              F(K,I) = F(K,I) + F1(K)
              FDOT(K,I) = FDOT(K,I) + F1DOT(K)
   25     CONTINUE
          PHI(I) = PHI(I) - BODY(J)/SQRT(A(1)**2 + A(2)**2 + A(3)**2)
*
          IF (J.EQ.KDUM) THEN
              J = J + 1
              GO TO 12
          END IF
   30 CONTINUE
   40 CONTINUE
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(I1,I2,1)
      END IF
*
*       Check case of new c.m. (KCASE = 1 with I1 = ICOMP, I2 = JCOMP).
      IF (KCASE.EQ.1) THEN
*       Form total c.m. force & first derivative.
          A1 = BODY(ICOMP)
          A2 = BODY(JCOMP)
          DO 60 K = 1,3
              F(K,NTOT) = (A1*F(K,ICOMP) + A2*F(K,JCOMP))/BODY(NTOT)
              FDOT(K,NTOT) = (A1*FDOT(K,ICOMP) + A2*FDOT(K,JCOMP))/
     &                                                        BODY(NTOT)
   60     CONTINUE
      END IF
*
*       Obtain new time-step for c.m. or single particles.
      IF (KCASE.EQ.1) THEN
          CALL STEPS(NTOT,NTOT,1)
      ELSE
          CALL STEPS(I1,I2,KCASE)
      END IF
*
      ISEND = -1
      RETURN
*
      END
