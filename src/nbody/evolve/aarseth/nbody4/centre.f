      SUBROUTINE CENTRE
*
*
*       Cluster centre diagnostics.
*       ---------------------------
*
      INCLUDE 'common4.h'
      REAL*8 R2,RHO
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
*
*       Ensure Lagrange radii are available if not done before.
      IF (KZ(7).EQ.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Determine central potential from five innermost members.
      PHI0 = 0.d0
      DO 10 I = 1,5
          IM = JLIST(I)
          PHI0 = PHI0 + PHI(IM)
   10 CONTINUE
      PHI0 = PHI0/5.d0
*
*       Form mean square velocity inside first Lagrangian radius (1 %).
      VC2 = 0.d0
      ZM = 0.d0
      ZM1 = 0.01d0*ZMASS
      I = 0
   15 I = I + 1
      IM = JLIST(I)
      VI2 = 0.d0
      DO 20 K = 1,3
          VI2 = VI2 + XDOT(K,IM)**2
   20 CONTINUE
      VC2 = VC2 + VI2
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZM1) GO TO 15
      VC2 = VC2/FLOAT(I)
*
*       Save radius of innermost Lagrange radius.
      RL1 = SQRT(R2(I))
*
*       Define central density (modulus 4*pi/3) inside 1 % mass percentile.
      RHOC = ZM/(R2(I)*SQRT(R2(I)))
*
*       Obtain alternative core radius from dynamical criterion.
      RC1 = SQRT(VC2/RHOC)
*
*       Evaluate core mass and membership inside RC1.
      ZMC1 = 0.d0
      NC1 = 0
      I = 0
   25 I = I + 1
      IM = JLIST(I)
      ZMC1 = ZMC1 + BODY(IM)
      NC1 = NC1 + 1
      IF (R2(I).LT.RC1**2) GO TO 25
*
*       Define energy of soft binaries (kT = 2/3 <m><v>**2; EB > -0.1*kT).
      ZKT = 0.6667d0*ZKIN/FLOAT(N - NPAIRS)
      SEMI0 = 0.5d0*BODYM**2/(0.1d0*ZKT)
      DTM = 0.02d0*SQRT(ETA/0.02d0)*SQRT(SEMI0**3/BODYM)
*
      NSB = 0
      EBS = 0.d0
      EBSM = 0.d0
*
*       Search for soft binaries using neighbourlist from GRAPE.
      DO 40 I = IFIRST,N
          IF (STEP(I).GT.DTM.OR.BODY(I).EQ.0.D0) GO TO 40
*       Define neighbour distance from density fitting expression.
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2
     &                                 + (X(3,I) - RDENS(3))**2
          RS2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
          CALL NBLIST(I,RS2,NNB)
          IF (NNB.EQ.0) GO TO 40
          JM = 0
          RIJM = 1.d0
          VIJM = 1.d0
          DO 30 L = 2,NNB+1
              J = ILIST(L)
              RIJ2 = 0.d0
              VIJ2 = 0.d0
              DO 28 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
   28         CONTINUE
              IF (RIJ2.LT.RIJM) THEN
                  JM = J
                  RIJM = RIJ2
                  VIJM = VIJ2
              END IF
   30     CONTINUE
          IF(JM.GT.0)THEN
              HI = 0.5d0*VIJM - (BODY(I) + BODY(JM))/SQRT(RIJM)
              EB = BODY(I)*BODY(JM)*HI/(BODY(I) + BODY(JM))
*       Include energies < -0.1*kT and only count each binary once.
              IF (EB.LT.-0.1*ZKT.AND.JM.GT.I) THEN
                  NSB = NSB + 1
                  EBS = EBS + EB
                  EBSM = MIN(EB,EBSM)
              END IF
          END IF
   40 CONTINUE
*
*       Find maximum KS binding energy and associated components.
      EBKM = 0.d0
      J1M = 1
      J2M = 1
      NEW = 0
      DO 50 IPAIR = 1,NPAIRS
          IF (BODY(N+IPAIR).GT.0.0) THEN
              J1 = 2*IPAIR - 1
              J2 = J1 + 1
              IF (LIST(2,J2).EQ.0) NEW = NEW + 1
              EBK = BODY(J1)*BODY(J2)*H(IPAIR)/BODY(N+IPAIR)
              IF (EBK.LT.EBKM) THEN
                  EBKM = EBK
                  J1M = J1
                  J2M = J2
              END IF
          END IF
   50 CONTINUE
*
      EBKM = MIN(EMERGE,EBKM)
      WRITE (6,60)
   60 FORMAT (/,6X,'T6    PHI0   RC1  NC1   MC1    RL1   NSB   EBS',
     &         '       EBSM     NKS  NEW    EBK      EBKM         NAME')
      I6 = TTOT*TSCALE
      WRITE (6,70)  I6, PHI0, RC1, NC1, ZMC1, RL1, NSB, EBS, EBSM,
     &              NPAIRS, NEW, EBIN+EMERGE, EBKM, NAME(J1M), NAME(J2M)
   70 FORMAT (' #8',I5,F8.3,F6.2,I5,2F7.3,I5,2F10.5,2I5,2F9.4,2I6)
*
      RETURN
*
      END
