      SUBROUTINE CSTAB4(ITERM)
*
*
*       Four-body chain stability test.
*       -------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*8  M,MB,MB1,MB2,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3),
     &        XCM2(3),VCM2(3)
      INTEGER  IJ(NMX)
*
*
*       Sort particle separations (K1 & K2 form closest pair).
      CALL R2SORT(IJ,R2)
*
*       Save indices with I1 & I2 as wide binary and I3 & I4 as outer body.
      I1 = IJ(3)
      I2 = IJ(4)
      I3 = IJ(1)
      I4 = IJ(2)
      MB = M(I1) + M(I2)
      MB2 = M(I3) + M(I4)
      MB1 = MB + MB2
*
*       Use perturbed triple stability instead for small middle distance.
      IF (1.0/RINV(2).LT.0.5*RSUM) THEN
          CALL CSTAB3(ITERM)
          GO TO 40
      END IF
*
*       Form orbital parameters with c.m. of I3 & I4 as third body.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      VREL34 = 0.0D0
      RDOT = 0.0D0
      RDOT3 = 0.0D0
      RB0 = 0.0
      RI2 = 0.0D0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          RB0 = RB0 + (X(J1) - X(J2))**2
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          XCM2(K) = (M(I3)*X(J3) + M(I4)*X(J4))/MB2
          VCM2(K) = (M(I3)*V(J3) + M(I4)*V(J4))/MB2
          RI2 = RI2 + (XCM2(K) - XCM(K))**2
          VREL21 = VREL21 + (VCM2(K) - VCM(K))**2
          VREL34 = VREL34 + (V(J3) - V(J4))**2
          RDOT3 = RDOT3 + (XCM2(K) - XCM(K))*(VCM2(K) - VCM(K))
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = XCM2(K)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = VCM2(K)
   10 CONTINUE
*
*       Evaluate orbital elements for inner and outer motion.
      RB = SQRT(R2(I1,I2))
      R3 = SQRT(RI2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      SEMI1 = 2.0/R3 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      ECC1 = SQRT((1.0D0 - R3/SEMI1)**2 + RDOT3**2/(SEMI1*MB1))
      RB2 = SQRT(R2(I3,I4))
      SEMI2 = 2.0/RB2 - VREL34/MB2
      SEMI2 = 1.0/SEMI2
*
*       Skip if innermost triple is not stable (including ECC > 1).
      RB0 = SQRT(RB0)
      SEMI0 = 2.0/RB0 - VREL2/MB
      SEMI0 = 1.0/SEMI0
      IF (SEMI0.LT.0.0.OR.ECC.GT.1.0) THEN
          ITERM = 1
          GO TO 40
      END IF
      ECC0 = SQRT((1.0 - RB0/SEMI0)**2 + RDOT**2/(SEMI0*MB))
      XFAC0 = (1.0 + M(I3)/MB)*(1.0 + ECC)/SQRT(1.0 - ECC)
      PCRIT0 = 2.8*XFAC0**0.4*SEMI0
      PMIN0 = SEMI*(1.0 - ECC)
      IF (PMIN0.LT.PCRIT0) THEN
          ITERM = 1
          GO TO 40
      END IF
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
      QL = MB/MB2
      Q1 = MAX(M(I2)/M(I1),M(I1)/M(I2))
      Q3 = QL**0.33333
      Q13 = Q1**0.33333
      AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
      EK = AR*SEMI*(1.0D0 + ECC)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Replace the EK criterion by the MA analytical stability formula.
      Q0 = MB2/MB
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q0)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 1.0 + Q0
      END IF
      FE = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
      PCRIT = PCRIT*(1.0 + 0.1*ABS(SEMI2/SEMI))
*
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XCM,VCM,ALPHA)
*
*       Include fudge factor for inclination effect.
      YFAC = 1.0 - 0.3*ALPHA/180.0
      PCRIT = YFAC*PCRIT
*
*       Check hierarchical stability condition for bound close pair.
      ITERM = 0
      IF (PMIN.GT.PCRIT.AND.SEMI.GT.0.0.AND.SEMI1.GT.0.0.AND.
     &    RB.GT.SEMI) THEN
          ITERM = -1
          WRITE (6,20)  ECC, ECC1, SEMI, SEMI1, PMIN, PCRIT, EK, ALPHA
   20     FORMAT (' QUAD HIARCH   E =',F6.3,'  E1 =',F6.3,
     &                     '  A =',1P,E8.1,'  A1 =',E8.1,'  PM =',E9.2,
     &                     '  PC =',E9.2,'  EK =',E9.2,'  IN =',0P,F6.1)
          RI = SQRT(CM(1)**2 + CM(2)**2 + CM(3)**2)
          EMAX = 0.0
          WRITE (81,30)  TIMEC, RI, NAMEC(I3), QL, Q1, ECC, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, ALPHA, EMAX
   30     FORMAT (F9.5,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
          CALL FLUSH(81)
      END IF
*
   40 RETURN
*
      END
