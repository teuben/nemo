      SUBROUTINE HISTAB(IPAIR,J,PMIN,RSTAB)
*
*
*       Hierarchical stability condition.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  XX(3,3),VV(3,3)
*
*
*       Define KS & c.m. indices.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
*
*       Set semi-major axis & eccentricity of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
*     ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*     ECC = SQRT(ECC2)
*
*       Form stability ratio (Eggleton & Kiseleva, Ap.J. 455, 640, 1995).
*     Q = BODY(I)/BODY(J)
*     Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
*     Q3 = Q**0.33333
*     Q13 = Q1**0.33333
*       Note sign error in third term of Eq. (2).
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
*       Save hierarchical stability separation in COMMON variable.
*     PCRIT = AR*SEMI*(1.0D0 + ECC)
*
*       Determine the outer eccentricity.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 5 K = 1,3
          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
          RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    5 CONTINUE
      RIJ = SQRT(RIJ2)
      A1 = 2.0/RIJ - VIJ2/(BODY(I) + BODY(J))
*
*       Exit on hyperbolic orbit.
      IF (A1.LT.0.0) THEN
          PMIN = 1.0
          RSTAB = 0.0
          GO TO 20
      END IF
*
      SEMI1 = 1.0/A1
      ECC2 = (1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*(BODY(I) + BODY(J)))
      ECC1 = SQRT(ECC2)
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Evaluate the stability formula without fudge factor.
      Q = BODY(J)/BODY(I)
      IF (ECC1.LT.1.0) THEN
*       ZFAC = (1.0 + ECC1)/(1.0 - ECC1)**0.182*(1.0 + Q)/(1.0 + ECC)**3
          XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 1.0 + Q
      END IF
*     PCRIT = 2.6*ZFAC**0.355*SEMI*(1.0 + ECC)
      FE  = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
*
*       Exit if stability value exceeds upper limit.
      IF (PCRIT.GT.1.5*PMIN.AND.J.LE.N) THEN
          RSTAB = PCRIT
          GO TO 20
      END IF
*
*       Choose the most dominant triple in case of two binaries.
      YFAC = 1.0
      IF (J.GT.N) THEN
          SEMI2 = -0.5D0*BODY(J)/H(J-N)
          SFAC = (1.0 + Q)**0.4*SEMI
          SFAC2 = (1.0 + BODY(I)/BODY(J))**0.4*SEMI2
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
          IF (SFAC2.GT.SFAC) THEN
              PCRIT = PCRIT*(SFAC2/SFAC)
          END IF
      END IF
*
*       Determine inclination for triple configuration.
      IF (J.LE.N) THEN
          DO 10 K = 1,3
              XX(K,1) = X(K,I1)
              XX(K,2) = X(K,I2)
              XX(K,3) = X(K,J)
              VV(K,1) = XDOT(K,I1)
              VV(K,2) = XDOT(K,I2)
              VV(K,3) = XDOT(K,J)
  10      CONTINUE
          CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
      ELSE
          ANGLE = 0.0
      END IF
*
*       Include fudge factor in the stability condition.
      YFAC = YFAC - 0.3*ANGLE/180.0
      PCRIT = YFAC*PCRIT
      RSTAB = PCRIT
*
   20 RETURN
*
      END
