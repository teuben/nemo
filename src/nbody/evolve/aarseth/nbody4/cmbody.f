      SUBROUTINE CMBODY(ENERGY,NSYS)
*
*
*       Formation of c.m. body by collision.
*       ------------------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (NMX=10,NMX4=4*NMX)
      COMMON/RCLOSE/  RIJ(4,4),RCOLL4,QPERI4,SIZE4(4),ECOLL4,IP(4)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
      REAL*8  CM(6),A0(3),A2(3)
      CHARACTER*8  WHICH1
*
*
*       Distinguish between chain and triple or quad case (ICH > 0 or = 0).
      IF (IPHASE.EQ.9) THEN
          ICH = 1
      ELSE
*       Activate collision indicator (otherwise done in CHTERM).
          ICH = 0
          IPHASE = 9
      END IF
*
*       Specify global indices of subsystem (membership: NSYS = 2 - 5).
      IF (NSYS.EQ.2) THEN
*
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          I = N + KSPAIR
          DT = 0.1*STEP(I)
          IF (DT.GT.2.4E-11) THEN
              TIME2 = TIME - TPREV
              CALL STEPK(DT,DTN)
              TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
              TIME = MIN(TBLOCK,TIME)
          ELSE
              TIME = MIN(T0(I) + STEP(I),TBLOCK)
          END IF
          TIME0 = TIME
*
*       Check for hierarchical configuration.
          I1 = 2*KSPAIR - 1
          I2 = I1 + 1
          JCL = 0
          NP1 = LIST(1,I1) + 1
          DO 5 L = 2,NP1
              J = LIST(L,I1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              A12 = 0.0
              A22 = 0.0
              A1A2 = 0.0
              DO 2 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) -XDOT(K,J))
                  K1 = K + 1
                  IF (K1.GT.3) K1 = 1
                  K2 = K1 + 1
                  IF (K2.GT.3) K2 = 1
                  A0(K) = (X(K1,I1)-X(K1,I2))*(XDOT(K2,I1)-XDOT(K2,I2))
     &                  - (X(K2,I1)-X(K2,I2))*(XDOT(K1,I1)-XDOT(K1,I2))
                  A2(K) = (X(K1,J) - X(K1,I))*(XDOT(K2,J) - XDOT(K2,I))
     &                  - (X(K2,J) - X(K2,I))*(XDOT(K1,J) - XDOT(K1,I))
                  A12 = A12 + A0(K)**2
                  A22 = A22 + A2(K)**2
                  A1A2 = A1A2 + A0(K)*A2(K)
    2         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  ECC1 = SQRT(ECC2)
                  RP1 = A1*(1.0 - ECC1)
                  SEMI = -0.5*BODY(I)/H(KSPAIR)
                  ECC = 1.0 - R(KSPAIR)/SEMI
                  RA = SEMI*(1.0 + ECC)
                  SR = RP1/RA
                  GA = 2.0*BODY(J)*(RA/RP1)**3/BODY(I)
*       Determine inclination (8 bins of 22.5 degrees).
                  FAC = A1A2/SQRT(A12*A22)
                  ANGLE = 360.0*ACOS(FAC)/TWOPI
                  WRITE (6,4)  KSPAIR, NAME(J), H(KSPAIR), ECC, SEMI,
     &                         A1, RP1, GA, ECC1, SR, ANGLE
    4             FORMAT (' HIERARCHY:   KS NMJ H E A0 A1 RP GA E1 SR ',
     &                    'IN',2I6,F7.0,F9.5,1P,4E9.1,0P,F6.2,F6.1,F7.1)
              END IF
*
*       Select closest single body inside 0.5*RMIN as KS component.
              IF (RIP.LT.0.5*RMIN.AND.J.LE.N) THEN
                  IF (JCL.GT.0) THEN
                      IF (RIP.GT.RIP0) GO TO 5
                      JCL = J
                      RIP0 = RIP
                  ELSE
                      JCL = J
                      RIP0 = RIP
                  END IF
              END IF
    5     CONTINUE
*
*       Search for evidence of recent regularization.
*         NAM1 = NAME(2*KSPAIR-1)
*         NAM2 = NAME(2*KSPAIR)
*         NNB = LISTD(1)
*         DO 7 K = 2,NNB+1
*             IF (LISTD(K).EQ.NAM1.OR.LISTD(K).EQ.NAM2) THEN
*                 WRITE (6,6)  NAM1, NAM2, LISTD(K), K
*   6             FORMAT (' KS REMNANT:    NAM LISTD K  ',3I6,I4)
*             END IF
*   7     CONTINUE
*
*       Predict body #JCL to current time in case of no collision.
          IF (JCL.GT.0) CALL XVPRED(JCL,-1)
*
*       Ensure orbit is at pericentre (perturbed hyperbolic case is OK).
          SEMI = -0.5*BODY(N+KSPAIR)/H(KSPAIR)
          IF (R(KSPAIR).GT.SEMI.AND.SEMI.GT.0.0) THEN
              CALL KSAPO(KSPAIR)
              CALL KSPERI(KSPAIR)
*       Restore quantized time to avoid small STEP (KSTERM needs T0 = TIME).
              TIME = TIME0
              T0(I1) = TIME
          END IF
*
*       Save collision distance and VINF before any common envelope stage.
          RCOLL = R(KSPAIR)
          VINF = 0.0
          IF (H(KSPAIR).GT.0.0) VINF = SQRT(2.0*H(KSPAIR))*VSTAR
          ECC = 1.0 - R(KSPAIR)/SEMI
*
*       Include diagnostics for every collision/coalescence event.
          IF (KZ(19).GE.3) THEN
              ZM1 = BODY(I1)*ZMBAR
              ZM2 = BODY(I2)*ZMBAR
              R1 = RADIUS(I1)*SU
              R2 = RADIUS(I2)*SU
              WRITE (86,9)  TPHYS, NAME(I1), NAME(I2), KSTAR(I1),
     &                      KSTAR(I2), ZM1, ZM2, R1, R2, RCOLL*SU, ECC
    9         FORMAT (' COLL:    TPH NAM K* M R* QP E ',
     &                           F8.1,2I6,2I4,2F6.2,3F7.1,F8.4)
              CALL FLUSH(86)
          END IF
*
*       Include special procedure for common envelope stage with mass loss.
          IF (KZ(19).GE.3) THEN
              K1 = KSTAR(I1)
              K2 = KSTAR(I2)
              ICASE = KTYPE(K1,K2)
              IF(ICASE.GT.100)THEN
*                 CALL EXPEL(I1,I2,ICASE)
                  IF (ICASE.LT.0) GO TO 100
*       Treat collision as before in case of CE without coalescence.
                  ICOMP = I1
              END IF
          END IF
*
*       Update body #JCL to current time for new KS with combined c.m.
          IF (JCL.GT.0) THEN
*             CALL XVPRED(JCL,-1)
              T0(JCL) = TIME
              CALL DTCHCK(TIME,STEP(JCL),DTK(40))
              DO 10 K = 1,3
                  X0DOT(K,JCL) = XDOT(K,JCL)
                  X0(K,JCL) = X(K,JCL)
   10         CONTINUE
          END IF
*
*       Check diagnostics of degenerate binary (skip case of velocity kick).
          IF(KZ(8).GT.3.AND.MAX(KSTAR(I1),KSTAR(I2)).GE.10)THEN
              IF(KSTAR(I1).LE.12)THEN
                  CALL DEGEN(KSPAIR,KSPAIR,5)
              END IF
          END IF
*
*       Save binding energy (BODY(I2) = 0 is OK).
          EB = BODY(I1)*BODY(I2)*H(KSPAIR)/BODY(I)
          WHICH1 = ' BINARY '
          IF (H(KSPAIR).GT.0.0) THEN
              WHICH1 = ' HYPERB '
              NHYP = NHYP + 1
          END IF
*
*       Terminate KS pair and set relevant indices for collision treatment.
          T0(I1) = TIME
          PHI1 = PHI(I)
          CALL DTCHCK(TIME,STEP(I1),DTK(40))
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
          I3 = 0
          ICOMP = I1
          DMIN2 = MIN(DMIN2,RCOLL)
      ELSE
*       Ignore case of three-body system here (JLIST(4) = 0).
          I1 = JLIST(1)
          I2 = JLIST(2)
          I3 = JLIST(3)
          I4 = JLIST(4)
          IQCOLL = 5
          VINF = 0.0
          PHI1 = -VC**2
          ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
          ECC = MAX(ECC,0.001D0)
          IF (EBS.GT.0) THEN
              HI = EBS*(BODY(I1) + BODY(I2))/(BODY(I1)*BODY(I2))
              VINF = SQRT(2.0*HI)*VSTAR
          END IF
*
*       Set new quantized time (note: restore if problems in CHTERM).
*         TIME = TBLOCK
*
*       Include special treatment for common envelope stage inside chain.
          IF (ICH.GT.0.AND.KZ(19).GE.3) THEN
              ZM1 = BODY(I1)*ZMBAR
              ZM2 = BODY(I2)*ZMBAR
              R1 = RADIUS(I1)*SU
              R2 = RADIUS(I2)*SU
              ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
              ECC = MAX(ECC,0.0D0)
              WRITE (86,9)  TPHYS, NAME(I1), NAME(I2), KSTAR(I1),
     &                      KSTAR(I2), ZM1, ZM2, R1, R2, DMINC*SU, ECC
*
              K1 = KSTAR(I1)
              K2 = KSTAR(I2)
              ICASE = KTYPE(K1,K2)
              IF(ICASE.GT.100)THEN
                  IQCOLL = 6
*                 CALL EXPEL2(I1,I2,ICASE)
*       Decide between chain restart, coalescence or collision.
                  IF (ICASE.GT.0) THEN
*       Adopt negative membership and reverse NSYS to denote chain restart.
                      NCH = -NSYS
                      NSYS = -NSYS
                      GO TO 100
                  END IF
*       Check for coalescence (one body of zero mass).
                  IF (BODY(I1).EQ.0.0D0.OR.BODY(I2).EQ.0.0D0) THEN
*       Reduce the membership (< 0 for SETSYS) and remove ghost from chain.
                      NCH = -(NSYS - 1)
                      JLIST(1) = I1
                      IF (BODY(I1).EQ.0.0D0) JLIST(1) = I2
                      DO 11 L = 2,NSYS-1
                          JLIST(L) = JLIST(L+1)
   11                 CONTINUE
                      GO TO 100
                  END IF
*       Treat collision as before in case of CE without coalescence.
                  ICOMP = I1
                  JLIST(1) = I1
                  JLIST(2) = I2
              END IF
          END IF
      END IF
*
*       Obtain mass loss and evolution epoch of composite star.
      DM = 0.0D0
      IF (KZ(19).GE.3) THEN
          CALL MIX(I1,I2,DM)
          ICOMP = I1
*       Note possible switching of I1 and I2 (cf. JLIST).
      END IF
*
*       Define global c.m. coordinates & velocities from body #I1 & I2.
      ZM = BODY(I1) + BODY(I2)
      DO 12 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
   12 CONTINUE
*
*	Ensure T0 = TIME for correct potential energy correction in FCORR.
      IF (ICH.GT.0) THEN
          DO 14 L = 1,NCH
              J = JLIST(L)
              T0(J) = TIME
              CALL DTCHCK(TIME,STEP(J),DTK(40))
   14     CONTINUE
      END IF
*
*       Evaluate potential energy with respect to colliding bodies.
      IF (NSYS.EQ.2) THEN
*       Form list of close neighbours and copy members to array JPERT.
          RS2 = (0.1*RC)**2
          CALL NBLIST(I1,RS2,NNB)
          DO 15 L = 1,NNB
              JPERT(L) = ILIST(L+1)
   15     CONTINUE
          JLIST(1) = I1
          JLIST(2) = I2
*       Replace second old KS component temporarily by arbitrary body.
          JPERT(1) = N
          IF (I2.EQ.N) JPERT(1) = N + 1
          CALL NBPOT(2,NNB,POT1)
      ELSE
*       Obtain differential effect on #I1 & #I2 due to other members.
          DO 16 L = 3,NCH
              JPERT(L-2) = JLIST(L)
   16     CONTINUE
          NP = NCH - 2
          CALL NBPOT(2,NP,POT1)
      END IF
*
*       Create new body from c.m. and initialize zero mass ghost in #I2.
      BODY(I1) = ZM
      BODY(I2) = 0.0D0
      SPIN(I1) = (SPIN(I1) + SPIN(I2))*(1.0 - DM/ZM)
      PHI(I1) = PHI1
      T0(I2) = TADJ + DTADJ 
*     STEP(I2) = 1.0D+06
      DTMAX = DTK(1)
      CALL DTCHCK(TIME,DTMAX,DTK(40))
      STEP(I2) = DTMAX
      RI = SQRT(X(1,I2)**2 + X(2,I2)**2 + X(3,I2)**2)
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
*
      DO 20 K = 1,3
          X(K,I1) = CM(K)
          X0(K,I1) = CM(K)
          XDOT(K,I1) = CM(K+3)
          X0DOT(K,I1) = CM(K+3)
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I2) = 1000.0*RSCALE*X(K,I2)/RI
          X(K,I2) = X0(K,I2)
          X0DOT(K,I2) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I2)/VI
          XDOT(K,I2) = X0DOT(K,I2)
          F(K,I2) = 0.0D0
          FDOT(K,I2) = 0.0D0
          D2(K,I2) = 0.0D0
          D3(K,I2) = 0.0D0
   20 CONTINUE
*
*       Update variables on GRAPE (needed for NBLIST if DM > 0).
      CALL GPSEND
*
*       Refresh index of dominant body in case of switch in routine MIX.
      JLIST(1) = I1
*       Obtain potential energy w.r.t. new c.m. and apply tidal correction.
      IF (NSYS.EQ.2) THEN
          CALL NBPOT(1,NNB,POT2)
      ELSE
          CALL NBPOT(1,NP,POT2)
      END IF
      DP = POT2 - POT1
      ECOLL = ECOLL + DP
*
*       Remove the ghost particle from perturber lists containing #I1.
      JPERT(1) = I2
      JLIST(1) = I2
      CALL NBREM(I1,1,NPAIRS)
*
*       Include correction procedure in case of mass loss (cf routine MIX).
      IF (KZ(19).GE.3.AND.DM.GT.0.0D0) THEN
*
*       Form neighbour list (depending on mass loss; IPHASE may be < 0).
          FACM = MIN(1.0 + DM*ZMBAR,3.0D0)
          RS2 = (FACM*RSCALE)**2/FLOAT(N)**0.66667
          CALL NBLIST(I1,RS2,NNB)
          ILIST(1) = NNB
*
*       Reduce mass of composite body and update total mass (check SN mass).
          BODY(I1) = ZM - DM
          BODY(I1) = MAX(BODY(I1),0.0D0)
          ZMASS = ZMASS - DM
*
*       Perform total force & energy corrections (new polynomial set later).
          KW = KSTAR(I1)
          CALL FCORR(I1,DM,KW)
*
*       Initialize new polynomials of neighbours & #I for DM > 0.1 DMSUN.
          IF (DM*ZMBAR.GT.0.1) THEN
*
*       Include body #I at the end (counting from location #2; not KS case).
              NNB2 = NNB + 2
              ILIST(NNB2) = I1
              IF (NSYS.EQ.2) NNB2 = NNB2 - 1
*       Specify IPHASE = -3 to preserve ILIST and skip GPSEND in FPOLYI.
              IPHASE = -3
              ISEND = -1
*
*       Obtain new F & FDOT and time-steps.
              DO 30 L = 2,NNB2
                  J = ILIST(L)
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
*                 CALL FPOLY1(J,J,0)
                  CALL FPOLYI(J)
   30         CONTINUE
          END IF
          ISEND = -1
          TPREV = TIME - STEPX
      END IF
*
*       Decide appropriate path for each case.
      IF (NSYS.EQ.2) GO TO 40
      IF (NSYS.EQ.3) GO TO 45
*
*       Switch KS components if body #I3 & I4 is closer than #I1 & I3.
      IF (JLIST(6).LT.0) THEN
          I4 = I1
          I1S = I1
          I1 = JLIST(4)
          JLIST(4) = I1S
      END IF
      ICOMP = I4
*
*       Check KS case for new regularization with close hierarchical body.
   40 IF (NSYS.EQ.2) THEN
          IF (JCL.GT.0) THEN
              ICOMP = I1
              JCOMP = JCL
              CALL KSREG
              GO TO 80
          END IF
      END IF
*
*       Initialize force polynomial for new single or third body (ICOMP).
*     CALL FPOLY1(ICOMP,ICOMP,0)
      ISEND = -1
      IPHASE = 0
      CALL FPOLYI(ICOMP)
      IF (NSYS.EQ.2) GO TO 80
*
*       Add kinetic energy from last body and check DMIN in TRIPLE or QUAD.
   45 IF (ICH.GT.0) THEN
          WHICH1 = '  CHAIN '
          RCOLL = DMINC
*       Specify new membership (< 0 for SETSYS) and remove ghost from chain.
          NCH = -(NSYS - 1)
          JLIST(1) = I1
          DO 50 L = 2,NSYS-1
              JLIST(L) = JLIST(L+1)
   50     CONTINUE
*       Copy well defined binding energy and skip explicit evaluation.
          EB = EBS
          CHCOLL = CHCOLL + EB
          GO TO 80
      ELSE IF (NSYS.EQ.3) THEN
          WHICH1 = ' TRIPLE '
          DMIN3 = MIN(DMIN3,RCOLL4)
          RCOLL = RCOLL4
      ELSE
          WHICH1 = '   QUAD '
          DMIN4 = MIN(DMIN4,RCOLL4)
          RCOLL = RCOLL4
      END IF
*
*       Obtain binding energy of the subsystem (ignore ghost).
      JLIST(1) = I1
      ZKE = 0.0D0
      POTS = 0.0D0
      DO 60 L = 1,NSYS
          I = JLIST(L)
          ZKE = ZKE + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                        XDOT(3,I)**2)
   60 CONTINUE
*
      DO 70 L = 1,NSYS-1
          DO 65 LL = L+1,NSYS
              I = JLIST(L)
              J = JLIST(LL)
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              POTS = POTS + BODY(I)*BODY(J)/SQRT(RIJ2)
   65     CONTINUE
   70 CONTINUE
*
*       Form net energy correction for triple or quad case.
      EB = ENERGY - (0.5D0*ZKE - POTS)
*
*       Set global components for new KS regularization (ICOMP < JCOMP).
      ICOMP = MIN(I1,I3)
      JCOMP = MAX(I1,I3)
*
*       Initialize new KS pair (set ISEND < 0 for FPOLYI).
      ISEND = -1
      CALL KSREG
*
*       Update energy loss & collision counters.
   80 ECOLL = ECOLL + EB
      E(10) = E(10) + EB + DP
*       Distinguish case of contact binary (i.e. coalescence).
      IF (IQCOLL.EQ.3) THEN
          NPOP(8) = NPOP(8) + 1
          NCOAL = NCOAL + 1
          WRITE (6,85)  IQCOLL, NAME1, NAME2, ZM, RCOLL, EB, DP, ECC
   85     FORMAT (/,' BINARY COAL    IQCOLL =',I3,'  NAME =',2I6,
     &             '  M =',F7.4,'  RCOLL =',1P,E8.1,' EB =',E9.1,
     &             '  DP =',E9.1,'  E =',0P,F8.4)
          GO TO 95
      END IF
*
      NPOP(8) = NPOP(8) + 1
      NCOLL = NCOLL + 1
*
      WRITE (6,90)  WHICH1, NSYS, NAME1, NAME2, ZM, RCOLL, EB, VINF,
     &              ECC, DP
   90 FORMAT (/,A8,'COLLISION    NSYS =',I3,'  NAME =',2I6,
     &             '  M =',F7.4,'  RCOLL =',1P,E8.1,'  EB =',E9.1,
     &             '  VINF =',0P,F5.1,'  ECC =',F9.5,'  DP =',1P,E9.1)
*
*       Specify IPHASE < 0 for new sorting and ISEND < 0 for new GPSEND.
   95 IPHASE = -1
      ISEND = -1
*
*       Reduce NSUB for chain (temporary increase by CHINIT before CHTERM).
  100 IF (ICH.GT.0) THEN
          NSUB = NSUB - 1
      END IF
      TTOT = TIME + TOFF
*
      RETURN
*
      END
