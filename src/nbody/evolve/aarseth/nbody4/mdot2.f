      SUBROUTINE MDOT2
*
*
*       Mass loss from evolving stars (Chernoff-Weinberg).
*       --------------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  M1,TEV1,DLOGM,LUM
*
*
*       Define current time (unit of million years) and update counter.
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
*
*       Find global index of next star to be checked (NB! reset KS & IKS).
    1 KS = 0
      IKS = 0
      ITRY = 0
      DO 10 J = 1,NTOT
          IF (TEV(J).LE.TIME) THEN
              I = J
              GO TO 20
          END IF
   10 CONTINUE
*
*       Determine new evolution time (current TMDOT may be escaped star).
      GO TO 70
*
*       Skip possible c.m. of compact subsystem or zero mass ghost.
   20 IF (NSUB.GT.0) THEN
          IF (BODY(I).GT.BODY1.OR.NAME(I).EQ.0.OR.BODY(I).LE.0.0) THEN
              TEV(I) = TEV(I) + 0.01d0*TCR
              GO TO 1
          END IF
      END IF
*
*       Terminate merger if #I is first or the second comp is a c.m. body.
      IF (I.LT.IFIRST) THEN
          KSPAIR = KVEC(I)
          IF ((NAME(N+KSPAIR).LT.0.AND.I.EQ.2*KSPAIR - 1).OR.
     &        (NAME(I).GT.NZERO.AND.I.EQ.2*KSPAIR)) THEN
              IPHASE = 7
              CALL RESET
              GO TO 70
          END IF
      END IF
*
*       Extend time-scale if body #I is a ghost member (merger or chain).
      IF (BODY(I).LE.0.0D0) THEN
          TEV(I) = TEV(I) + 0.01d0*TCR
          GO TO 1
      END IF
*
*       Specify Chernoff-Weinberg mass loss and set large evolution time.
      NMDOT = NMDOT + 1
      ZM = BODY(I)*ZMBAR
      IF (ZM.LT.4.7) THEN
          DM = ZM - (0.58d0 + 0.22d0*(ZM - 1.d0))
*         DM = MAX(DM,0.2d0*ZM)
          KW = 10
          NWD = NWD + 1
          RM = 0.01d0
          LUM = 2.d0
      ELSE IF (ZM.GT.8.0) THEN
          DM = ZM - 1.4d0
          KW = 13
          NSN = NSN + 1
          RM = 1.0d-05
          LUM = 0.03d0
      ELSE
          DM = ZM
          KW = 0
          NRS = NRS + 1
          RM = 0.d0
          LUM = 0.d0
      END IF
      TEV(I) = 1.0d+10
      KW0 = KSTAR(I)
      KSTAR(I) = KW
      RNEW = RM/SU
      IF (RADIUS(I).GT.0.0) RADIUS(I) = RM/SU
      IF (ZLMSTY(I).GT.0.0) ZLMSTY(I) = LUM
*
*       Set spin to zero. 
      SPIN(I) = 0.D0
*
*       Set mass loss in scaled units and form relative mass ratio.
      DM = DM/ZMBAR
      DMR = ABS(DM/BODY(I))
      IF (I.GT.N) THEN
          WRITE (6,29) I, IFIRST, KW, KW0, DMR, DM*ZMBAR
   29     FORMAT (' DANGER!    MDOT:    I I* KW K* DMR DMS ',
     &                         I5,3I4,F7.3,1P,E10.2)
          TEV(I) = 1.0d+10
          GO TO 70
      END IF
*
*       Define KS index & Roche indicator and update core mass (chaos only).
      IF (I.LT.IFIRST) THEN
          KSPAIR = KVEC(I)
          IF (DMR.LT.0.01) THEN
              IF (NAME(N+KSPAIR).GT.0) ITRY = 2
          END IF
      END IF
*
*       Include special procedures for KS components.
      IF (I.LT.IFIRST.AND.DMR.GT.0.01) THEN
*       Distinguish between KS pair and merger configuration.
          SEMI = -0.5d0*BODY(N+KSPAIR)/H(KSPAIR)
          IF (NAME(N+KSPAIR).GT.0) THEN
*       Set random phase for neutron star formation with negative KS index.
              IF (KW.EQ.13) THEN
                  JPAIR = -KSPAIR
                  CALL KSAPO(JPAIR)
              END IF
*       Terminate for large mass loss or soft binary and re-determine index.
              IF ((DMR.GT.0.2.AND.R(KSPAIR).GT.RMIN).OR.
     &            H(KSPAIR) + DM/SEMI.GT.-ECLOSE.OR.KW.GE.13) THEN
                  I = I + 2*(NPAIRS - KSPAIR)
*       Predict current KS variables and save at end of routine RESOLV.
                  CALL RESOLV(KSPAIR,3)
                  IPHASE = 2
                  JCOMP = 0
                  CALL KSTERM
                  KS = 1
              ELSE
*       Implement mass loss and expand KS orbit at constant eccentricity.
                  CALL HCORR(I,DM,RNEW)
                  ITRY = 1
              END IF
          ELSE
*       Adopt KS treatment for single outer component or terminate merger.
              IF (I.EQ.2*KSPAIR.AND.NAME(I).LE.NZERO.AND.
     &            H(KSPAIR) + DM/SEMI.LT.-ECLOSE.AND.KW.NE.13) THEN
                  CALL HCORR(I,DM,RNEW)
              ELSE
                  IPHASE = 7
                  CALL RESET
                  GO TO 70
              END IF
          END IF
      END IF
*
*       Perform neighbour force corrections if mass loss is significant.
      IF (DMR.GT.0.01) THEN
*
*       Include optional diagnostics for mass loss orbit (filename MDOT).
          IF (KZ(21).GT.2) THEN
              CALL MTRACE(I,DM)
          END IF
*
*       Update mass and previous evolution time TEV0.
          BODY(I) = BODY(I) - DM
          TEV0(I) = TEV(I)
*
*       Accumulate total mass loss (solar units) and reduce cluster mass.
          DMSUN = DM*ZMBAR
          ZMDOT = ZMDOT + DMSUN
          ZMASS = ZMASS - DM
*
*       Check optional diagnostics.
          IF (KZ(19).GT.5) THEN
              WRITE (6,40)  I, NAME(I), KW, KW0, BODY(I)*ZMBAR,
     &                      DMSUN, ZMDOT, TPHYS
   40         FORMAT (' MDOT:    I NM KW K* MS DMS ZMDOT T6 ',
     &                           4I5,F6.1,F7.2,2F8.1)
          END IF
*
*       Replace any KS components by corresponding c.m. for main procedures.
          IF (I.LT.IFIRST) THEN
              IKS = I
              I = N + KSPAIR
              I1 = 2*KSPAIR - 1
*       Predict coordinates & velocities of any unperturbed KS components.
              IF (LIST(1,I1).EQ.0) THEN
                  CALL RESOLV(KSPAIR,1)
              END IF
          END IF
*
*       Form neighbour list (depending on mass loss; IPHASE may be < 0).
          FACM = MIN(1.d0 + DMSUN,3.d0)
          RS2 = 0.5d0*(FACM*RSCALE)**2/FLOAT(N)**0.66667
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
*       Adopt neighbour distance from interparticle spacing or density fit.
          RS2 = MAX(RS2,H2)
   52     CALL NBLIST(I,RS2,NNB)
          IF (NNB.EQ.0) THEN
              RS2 = 1.2*RS2
              GO TO 52
          END IF
          ILIST(1) = NNB
*
*       Perform total force & energy corrections (delay dF if DMSUN < 0.1).
          CALL FCORR(I,DM,KW)
*
*       Initialize new polynomials of neighbours & #I for DMSUN > 0.1.
          IF (DMSUN.GT.0.1) THEN
*
*       Include body #I at the end (counting from location #2).
              NNB2 = NNB + 2
              ILIST(NNB2) = I
*       Specify IPHASE = -3 for using current members of ILIST in FPOLYI.
              IPHASE = -3
*
*       Obtain new F & FDOT and time-steps.
              DO 60 L = 2,NNB2
                  J = ILIST(L)
                  IF (L.EQ.NNB2) J = I
                  DO 55 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   55             CONTINUE
                  ISEND = -1
                  CALL FPOLYI(J)
   60         CONTINUE
          END IF
*       Reduce TPREV and specify IPHASE < 0 for new time-step sorting.
          TPREV = TIME - STEPX
          IPHASE = -3
      END IF
*
*       Place any zero-mass ghosts at large distance for escape removal.
      IF (KW.EQ.0) THEN
          IF (IKS.GT.0) THEN
              CALL KSTERM
              I = IFIRST
              IF (ABS(BODY(I+1)).LT.1.0D-10) I = I + 1
          END IF
          BODY(I) = 0.d0
          T0(I) = 1.0d+06
          DO 65 K = 1,3
              X0DOT(K,I) = 0.d0
              XDOT(K,I) = 0.d0
              F(K,I) = 0.d0
              FDOT(K,I) = 0.d0
              D2(K,I) = 0.d0
              D3(K,I) = 0.d0
   65     CONTINUE
*       Set large X0 & X to avoid perturber selection and become escaper.
          X0(1,I) = 10.d0*RTIDE
          X(1,I) = 10.d0*RTIDE
*       Initialize polynomial & time-step for the massive component.
          I2 = IFIRST + 2 - I
          CALL GPSEND
          CALL FPOLYI(I2)
          ISEND = -1
          GO TO 70
      END IF
*
      IF (IKS.GT.0) THEN
*       Restore index in case of KS component (used as c.m. above).
          I = IKS
*       Re-initialize KS polynomials for perturbed motion.
          IF (LIST(1,I1).GT.0) THEN
              CALL RESOLV(KSPAIR,1)
              CALL KSPOLY(KSPAIR,1)
          END IF
      END IF
*
*       See if former KS pair can be regularized again (ensure ISEND < 0).
      IF (KS.GT.0) THEN
          ICOMP = IFIRST
          JCOMP = IFIRST + 1
          RIJ2 = (X(1,ICOMP) - X(1,JCOMP))**2 +
     &           (X(2,ICOMP) - X(2,JCOMP))**2 +
     &           (X(3,ICOMP) - X(3,JCOMP))**2
          IF (RIJ2.LT.RMIN22) THEN
              ISEND = -1
*       Enforce a new neighbour list in routine FPOLYI.
              IPHASE = 1
              CALL KSREG
              IPHASE = -3
          END IF
      END IF
*
*       Determine the time for next stellar evolution check.
   70 TMDOT = 1.0d+10
      DO 80 J = 1,NTOT
          IF (TEV(J).LE.TMDOT) THEN
              TMDOT = TEV(J)
          END IF
   80 CONTINUE
*
*       Update maximum mass and turnoff mass (skip compact subsystems).
      IF (NSUB.EQ.0) THEN
          BODY1 = 0.d0
          DO 90 J = 1,N
              BODY1 = MAX(BODY1,BODY(J))
   90     CONTINUE
          TURN = BODY1*ZMBAR
      END IF
*
*       Check optional updating of TEV & TMDOT for all masses.
      IF (KZ(13).NE.0) THEN
          DO 100 I = 1,N
*             IF (BODY(I).GT.0.9*BODY1) THEN
                  M1 = BODY(I)*ZMBAR
                  if (kstar(i).eq.10.or.kstar(i).eq.13) then
                     tev(i) = 1.0d+10
                  else
                     TEV1 = TEV(I)
                     CALL INTEV(M1,TEV1,DLOGM)
                     TEV(I) = TEV1/TSTAR
                     IF (TEV(I).LE.TMDOT) THEN
                        TMDOT = TEV(I)
                     END IF
                  END IF
  100     CONTINUE
      END IF
*
      RETURN
*
      END
