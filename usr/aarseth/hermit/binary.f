      SUBROUTINE BINARY
*
*
*       Initial binary.
*       ---------------
*
      INCLUDE 'commonp.h'
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3)
*
*
*       Introduce binary components from relative motion.
          I = 2
          ZI = 1.0
*
*       Randomize perihelion, node & inclination.
          A1 = 0.0
          A2 = 0.25
          PI = TWOPI*A1
          OMEGA = TWOPI*A2
          ZI = 0.25*TWOPI*ZI
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PI)*COS(OMEGA) - SIN(PI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PI)*COS(OMEGA) - COS(PI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PI)*SIN(OMEGA) + SIN(PI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PI)*SIN(OMEGA) + COS(PI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PI)*SIN(ZI)
          QX(3) = COS(PI)*SIN(ZI) 
*
*       Specify component masses (copy BODY0 from IMF2 or use RATIO).
          I1 = 2*I - 1
          I2 = 2*I
          IF (KZ(20).GE.2) THEN
              BODY(I1) = BODY0(I1)
              BODY(I2) = BODY0(I2)
          ELSE IF (RATIO.EQ.1.0) THEN
              BODY(I2) = BODY(I1) 
          ELSE
              BODY(I1) = RATIO*BODY(I1)
              BODY(I2) = BODY(I1)*(1.0 - RATIO)/RATIO
          END IF
*
*       Choose random (thermalized) or fixed eccentricity.
          IF (ECC0.LT.0.0) THEN
              ECC2 = RAN2(IDUM1)
              ECC = SQRT(ECC2)
          ELSE
              ECC = ECC0
          END IF
*
*       Select semi-major axis from uniform distribution in log(A) or SEMI0.
          IF (RANGE.GT.0.0) THEN
              EXP = RAN2(IDUM1)*LOG10(RANGE)
              SEMI = SEMI0/10.0**EXP
          ELSE
              SEMI = SEMI0
          END IF
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMI*(1.0 + ECC)
          XORB(2) = 0.0
          VORB(1) = 0.0
          ZMBIN = BODY(I1) + BODY(I2)
          VORB(2) = SQRT(ZMBIN*(1.0D0 - ECC)/(SEMI*(1.0D0 + ECC)))
          EBIN0 = EBIN0 - 0.5*BODY(I1)*BODY(I2)/SEMI
*
*       Transform to relative variables.
          DO 40 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   40     CONTINUE
*
*       Set global variables for each component.
          DO 50 K = 1,3
              X(K,I1) = X(K,I1) + BODY(I2)*XREL(K)/ZMBIN
              X(K,I2) = X(K,I1) - XREL(K)
              XDOT(K,I1) = XDOT(K,I1) + BODY(I2)*VREL(K)/ZMBIN
              XDOT(K,I2) = XDOT(K,I1) - VREL(K)
   50     CONTINUE
   60 CONTINUE
*
      RETURN
*
      END
