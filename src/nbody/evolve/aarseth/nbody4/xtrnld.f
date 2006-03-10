      SUBROUTINE XTRNLD(I1,I2,KCASE)
*
*
*       External force & derivative.
*       ----------------------------
*
      INCLUDE 'common4.h'
      REAL*8  FI(3),FD(3)
*
*
*       See whether to include the galactic tidal force.
      IF (KZ(14).EQ.1) THEN
          IF (TIDAL(1).GT.0.0.AND.KCASE.EQ.1) THEN
*       Include tidal force & first derivative (I1 = I2 for single body).
              DO 10 I = I1,I2
                  F(1,I) = F(1,I) + TIDAL(1)*X(1,I) + TIDAL(4)*XDOT(2,I)
                  F(2,I) = F(2,I) - TIDAL(4)*XDOT(1,I)
                  F(3,I) = F(3,I) + TIDAL(3)*X(3,I)
                  FDOT(1,I) = FDOT(1,I) + TIDAL(1)*XDOT(1,I) +
     &                                    TIDAL(4)*F(2,I)
                  FDOT(2,I) = FDOT(2,I) - TIDAL(4)*F(1,I)
                  FDOT(3,I) = FDOT(3,I) + TIDAL(3)*XDOT(3,I)
   10         CONTINUE
          END IF
      ELSE IF(KZ(14).GT.1) THEN
          DO 11 K = 1,3
             FI(K) = 0.0
             FD(K) = 0.0
   11     CONTINUE
          CALL XTRNLF(I1,FI,FD)
          DO 12 K = 1,3
              IF(KCASE.EQ.10)THEN
                  F(K,I1) = F(K,I1) + 0.5*FI(K)
                  FDOT(K,I1) = FDOT(K,I1) + ONE6*FD(K)
              ELSE IF(KCASE.EQ.11) THEN
                  F(K,I1) = F(K,I1) - 0.5*FI(K)
                  FDOT(K,I1) = FDOT(K,I1) - ONE6*FD(K)
              ELSE
                  F(K,I1) = F(K,I1) + FI(K)
                  FDOT(K,I1) = FDOT(K,I1) + FD(K)
             END IF
   12     CONTINUE
      END IF
*
      RETURN
*
      END
