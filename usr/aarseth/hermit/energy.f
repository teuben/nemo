      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'commonp.h'
*
*
*       Calculate the potential energy.
      ZKIN = 0.0D0
      POT = 0.0D0
      I = 1
   20 JMIN = I + 1
      POTJ = 0.0D0
*
      DO 30 J = JMIN,NMASS
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          POTJ = POTJ + BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3)
   30 CONTINUE
*
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.NMASS) GO TO 20
*
*       Sum the kinetic energy and solar interaction.
      DO 40 I = 1,NMASS
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
      ZKIN = 0.5D0*ZKIN
*
*       Total energy = ZKIN - POT.
*
      RETURN
*
      END
