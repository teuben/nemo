      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common1.h'
*
*
*       Calculate the potential energy.
      ZKIN = 0.0
      POT = 0.0
      I = 1
   10 POTJ = 0.0
      DO 15 J = I+1,N
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          POTJ = POTJ + BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3 + EPS2)
   15 CONTINUE
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.N) GO TO 10
*
*       Obtain the kinetic energy.
      DO 20 I = 1,N
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   20 CONTINUE
*
      ZKIN = 0.5*ZKIN
*       Total energy = ZKIN - POT.
*
      RETURN
*
      END
