      SUBROUTINE PHIDOT(I,DPHI)
*
*
*       Change in orbital potential.
*       ----------------------------
*
      INCLUDE 'common4.h'
*
*
*       Obtain integrated contribution from v*F.
      DPHI = 0.0
      DT = TIME - T0(I)
      DO 5 K = 1,3
*         DPHI = DPHI - (X0DOT(K,I) + 2.0*F(K,I)*DT)*
*    &                  (2.0*F(K,I) + 3.0*FDOT(K,I)*DT)
          DPHI = DPHI - 2.0*X0DOT(K,I)*F(K,I)
    5 CONTINUE
*
*       Subtract external tidal effect (already contained in F).
      IF (TIDAL(1).GT.0.0D0) THEN
          DO 10 K = 1,3
              DPHI = DPHI + TIDAL(K)*X0(K,I)*X0DOT(K,I)
   10     CONTINUE
      END IF
*
*       Form change in potential due to orbital motion since T0.
      DPHI = DPHI*DT
*
      RETURN
*
      END
