      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common2.h'
*
*
*       Calculate the potential & virial energy.
      ZKIN = 0.0
      POT = 0.0
      VIR = 0.0
      I = 1
   10 POTJ = 0.0
      VIRJ = 0.0
      DO 15 J = I+1,N
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DIJ2 = 1.0/(RIJ2 + EPS2)
          POTIJ = BODY(J)*SQRT(DIJ2)
          POTJ = POTJ + POTIJ
          VIRJ = VIRJ + POTIJ*RIJ2*DIJ2
   15 CONTINUE
      POT = POT + BODY(I)*POTJ
      VIR = VIR + BODY(I)*VIRJ
      I = I + 1
      IF (I.LT.N) GO TO 10
*
*       Obtain the kinetic energy.
      DO 20 I = 1,N
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   20 CONTINUE
      ZKIN = 0.5*ZKIN
*
*       Include optional external potential & virial (ETIDE & VIR > 0).
      ETIDE = 0.0
      IF (KZ(15).GT.0) THEN
          DO 30 I = 1,N
              RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
*       Distinguish between Plummer and logarithmic potential.
              IF (KZ(15).EQ.1) THEN
                  DI2 = RI2 + XTPAR(2)
                  ETIDE = ETIDE + BODY(I)*XTPAR(1)/SQRT(DI2)
                  VIR = VIR + BODY(I)*XTPAR(1)*RI2/(DI2*SQRT(DI2))
              ELSE
                  RI = SQRT(RI2)
                  ETIDE = ETIDE - CGAS*BODY(I)*LOG(RI/RGAS)
                  VIR = VIR + CGAS*BODY(I)
              END IF
   30     CONTINUE
      END IF
*
*       Total energy = ZKIN - POT - ETIDE.
*
      RETURN
*
      END
