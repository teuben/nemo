      SUBROUTINE STEPI(I)
*
*
*       Irregular time-step.
*       --------------------
*
      INCLUDE 'common2.h'
      REAL*4  F1DOT(4),F2DOT(4),F3DOT(4)
*
*
*       Set time-step differences.
      DT = TIME - T1(I)
      DT1 = TIME - T2(I)
      SI = DT + DT1
*
*       Form Taylor series force derivatives.
      DO 195 K = 1,3
          F1DOT(K) = D2(K,I)*DT + D1(K,I)
          F2DOT(K) = D3(K,I)*SI + D2(K,I)
          F3DOT(K) = D3(K,I)
  195 CONTINUE
*
*       Obtain new irregular integration step using composite expression.
      FI2 = FI(1,I)**2 + FI(2,I)**2 + FI(3,I)**2
      F1DOT(4) = F1DOT(1)**2 + F1DOT(2)**2 + F1DOT(3)**2
      F2DOT(4) = 4.0*(F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2)
      F3DOT(4) = 36.0*(F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2)
*
*       Prevent overflow for small softening (extra SQRT for large F3DOT).
      IF (F3DOT(4).LT.1.0E+20) THEN
          A1 = (SQRT(FI2*F2DOT(4)) + F1DOT(4))/
     &                              (SQRT(F1DOT(4)*F3DOT(4)) + F2DOT(4))
      ELSE
          A1 = (SQRT(FI2*F2DOT(4)) + F1DOT(4))/
     &                        (SQRT(F1DOT(4))*SQRT(F3DOT(4)) + F2DOT(4))
      END IF
*
*       Restrict increase of irregular time-step by stability factor 1.2.
      STEP(I) = MIN(SQRT(ETAI*A1),1.2*STEP(I),STEPR(I))
*
      RETURN
*
      END
