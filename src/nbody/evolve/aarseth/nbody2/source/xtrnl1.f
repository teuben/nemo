      SUBROUTINE XTRNL1(I,KCASE,FREG)
*
*
*       External force & derivative for Plummer model.
*       ----------------------------------------------
*
      INCLUDE 'common2.h'
      REAL*8  FREG(3)
*
*
*       Plummer potential: PHI(R) = -G*M/SQRT(R**2 + A**2).
*       Force component:   F(X) = -G*M*X/(R**2 + A**2)**1.5.
*
*       Form the inverse softened distance and M/R**3.
      RI2 = 1.0/(X(1,I)**2 + X(2,I)**2 + X(3,I)**2 + XTPAR(2))
      RI = SQRT(RI2)
      XF = XTPAR(1)*RI2*RI
*
*       Check call index (KCASE = 1: routine REGINT; = 2: FPOLY1).
      IF (KCASE.EQ.1) THEN
*       Add contribution to the regular force.
          DO 1 K = 1,3
              FREG(K) = FREG(K) - XF*X(K,I)
    1     CONTINUE
      ELSE
*       Obtain regular force & first derivative for new polynomial.
          RRDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
          XFDOT = -3.0*RRDOT*RI2
          DO 5 K = 1,3
              FR(K,I) = FR(K,I) - XF*X(K,I)
              D1R(K,I) = D1R(K,I) - XF*(XDOT(K,I) + XFDOT*X(K,I))
    5     CONTINUE 
      END IF
*
      RETURN
*
      END
