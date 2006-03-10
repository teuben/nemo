      SUBROUTINE XTRNLV(I1,I2)
*
*
*       External potential.
*       -------------------
*
      INCLUDE 'common4.h'
*
*
*       See whether to include the galactic tidal force.
      ET = 0.0D0
      IF (KZ(14).LE.2) THEN
          DO 10 I = I1,I2
              ET = ET - 0.5D0*BODY(I)*(TIDAL(1)*X(1,I)**2 +
     &                                 TIDAL(3)*X(3,I)**2)
   10     CONTINUE
      END IF
*
*       Place sum in ETIDE and single particle contribution in HT.
      IF (I2.GT.I1.AND.KZ(14).LE.2) THEN
          ETIDE = ET
      ELSE
          HT = ET
      END IF
*
      RETURN
*
      END
