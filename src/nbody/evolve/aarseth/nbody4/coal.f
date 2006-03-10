      SUBROUTINE COAL(IPAIR,KW1,KW2,MASS)
*
*
*       Coalescence of Roche/CE binary.
*       -------------------------------
*
      INCLUDE 'common4.h'
      CHARACTER*8  WHICH1
      REAL*8  CM(6)
      REAL*8  MASS(2)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Distinguish between KS and chain regularization.
      IF (IPAIR.GT.0) THEN
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          I = N + IPAIR
          DT = 0.1d0*STEP(I)
          IF (DT.GT.2.4E-11) THEN
              TIME2 = TIME - TPREV
              CALL STEPK(DT,DTN)
              TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
              TIME = MIN(TBLOCK,TIME)
          ELSE
              TIME = MIN(T0(I) + STEP(I),TBLOCK)
          END IF
*
*       Set zero energy (EB correction done in routines EXPEL & CHCORR).
          EB = 0.d0
          RCOLL = R(IPAIR)
          DMIN2 = MIN(DMIN2,RCOLL)
          VINF = 0.0
*
*       Define indicator for different cases, including hyperbolic KS.
          IF (KSTAR(I).LE.10) THEN
              WHICH1 = ' BINARY '
              IQCOLL = 0
          ELSE
              WHICH1 = '  ROCHE '
*       Save energy for correction (Roche COAL but not via CMBODY & EXPEL).
              IF (IQCOLL.EQ.0) THEN
                  EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
              END IF
              IQCOLL = 3
          END IF
          IF (H(IPAIR).GT.0.0) THEN
              WHICH1 = ' HYPERB '
              NHYP = NHYP + 1
              IQCOLL = -1
              VINF = SQRT(2.0*H(IPAIR))*VSTAR
          END IF
*
*       Terminate KS pair and set relevant indices for collision treatment.
          IPHASE = 9
          KSPAIR = IPAIR
          T0(2*IPAIR-1) = TIME
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
      ELSE
*       Copy dominant indices, two-body separation and binding energy.
          I1 = JLIST(1)
          I2 = JLIST(2)
          RCOLL = DMINC
          EB = 0.d0
          IQCOLL = 5
          WHICH1 = '  CHAIN '
*       Note that new chain TIME already quantized in routine CHTERM.
      END IF
*
*       Define global c.m. coordinates & velocities from body #I1 & I2.
      ICOMP = I1
      ZM = BODY(I1) + BODY(I2)
      DO 5 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
    5 CONTINUE
*
*       Form central distance (scaled by RC), central distance and period.
      RI2 = 0.d0
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      DO 10 K = 1,3
          RI2 = RI2 + (X(K,I1) - RDENS(K))**2
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
   10 CONTINUE
      RI = SQRT(RI2)/RC
      RIJ = SQRT(RIJ2)
      SEMI = 2.d0/RIJ - VIJ2/ZM
      SEMI = 1.d0/SEMI
      TK = DAYS*SEMI*SQRT(ABS(SEMI)/ZM)
*
*       Form perturber list (less #I2) and copy members to array JPERT.
      RS2 = (0.1d0*RC)**2
      CALL NBLIST(I1,RS2,NNB)
      L2 = 1
      JMIN = N + 1
      RIJ2 = 1.0d+10
      DO 15 L = 1,NNB
          JPERT(L) = ILIST(L+1)
          IF (JPERT(L).EQ.I2) THEN
              L2 = L
              JPERT(L) = N
              IF (I2.EQ.N) JPERT(L) = N - 1
          ELSE
              JJ = JPERT(L)
              RI2 = 0.d0
              DO 16 K = 1,3
                 RI2 = RI2 + (CM(K) - X(K,JJ))**2
 16           CONTINUE
              IF(RI2.LT.RIJ2)THEN
                  JMIN = JJ
                  RIJ2 = RI2
              ENDIF
          END IF
   15 CONTINUE
      IF(IPAIR.LE.0) JMIN = N + 1
*
*       Evaluate potential energy with respect to colliding bodies.
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NNB,POT1)
*
*       Specify new mass from sum and initialize zero mass ghost in #I2.
      ZMNEW = (MASS(1) + MASS(2))/ZMBAR
      DM = ZM - ZMNEW
      IF (DM.LT.1.0D-10) DM = 0.d0
*       Delay inclusion of any mass loss until after energy correction.
      M1 = BODY(I1)
      M2 = BODY(I2)
      BODY(I1) = ZM
      BODY(I2) = 0.d0
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
      IF(BODY0(I1).LT.BODY0(I2))THEN
         BODY0(I1) = BODY0(I2)
         EPOCH(I1) = EPOCH(I2)
         TEV(I1) = TEV(I2)
         SPIN(I1) = SPIN(I2)
         RADIUS(I1) = RADIUS(I2)
         NAME2 = NAME(I2)
         NAME(I2) = NAME(I1)
         NAME(I1) = NAME2
      ENDIF
      T0(I1) = TIME
      T0(I2) = TADJ + DTADJ 
      CALL DTCHCK(TIME,STEP(I2),DTK(40))
*
* Start the new star from the current time unless it has come from
* ROCHE and has TEV0 > TIME.
*
      TEV(I1) = MAX(TIME,TEV0(I1))
      TEV0(I1) = TEV(I1)
      TEV(I2) = 1.0d+10
      RI = SQRT(X(1,I2)**2 + X(2,I2)**2 + X(3,I2)**2)
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
*
*	Set T0 = TIME for any other chain members (cf. CALL POTI in FCORR).
      IF (IPAIR.LT.0) THEN
         DO 18 L = 1,NCH
            J = JLIST(L)
            IF (J.NE.I1.AND.J.NE.I2) THEN
               T0(J) = TIME
            END IF
   18    CONTINUE
      END IF
*
*       Check that a mass-less primary has type 15 for kick velocity.
      IF(ZMNEW*ZMBAR.LT.0.001.AND.KW1.NE.15)THEN
         WRITE(6,*)' ERROR COAL: mass1 = 0.0 and kw1 is not equal 15'
         WRITE(6,*)' I KW ',I1,KW1
         CALL gpfree
         STOP
      END IF
*
      DO 20 K = 1,3
          X(K,I1) = CM(K)
          X0(K,I1) = CM(K)
          XDOT(K,I1) = CM(K+3)
          X0DOT(K,I1) = CM(K+3)
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I2) = 1000.d0*RSCALE*X(K,I2)/RI
          X(K,I2) = X0(K,I2)
          X0DOT(K,I2) = SQRT(0.004d0*ZMASS/RSCALE)*XDOT(K,I2)/VI
          XDOT(K,I2) = X0DOT(K,I2)
          F(K,I2) = 0.d0
          FDOT(K,I2) = 0.d0
          D2(K,I2) = 0.d0
          D3(K,I2) = 0.d0
   20 CONTINUE
*
*       Obtain potential energy w.r.t. new c.m. and apply tidal correction.
      CALL NBPOT(1,NNB,POT2)
      ECOLL = ECOLL + (POT2 - POT1)
      DP = POT2 - POT1
      JPERT(L2) = I2
*
*       Remove the ghost particle from perturber lists containing #I1.
      JLIST(1) = I2
      CALL NBREM(I1,1,NPAIRS)
      JLIST(1) = I1
*
*       Include correction procedure in case of mass loss (cf routine MIX).
      IF (KZ(19).GE.3.AND.DM.GT.0.0) THEN
*
*       Reduce mass of composite body and update total mass (check SN mass).
          BODY(I1) = ZMNEW
          BODY(I1) = MAX(BODY(I1),0.d0)
          IF (ABS(BODY(I1)).LT.1.0d-10) TEV(I1) = 1.0d+10
          ZMASS = ZMASS - DM
*
*       Update variables on GRAPE.
          CALL GPSEND
*
*       Form neighbour list (depending on mass loss; IPHASE may be < 0).
          FACM = MIN(1.d0 + DM*ZMBAR,3.d0)
          RS2 = (FACM*RSCALE)**2/FLOAT(N)**0.66667
          CALL NBLIST(I1,RS2,NNB)
          ILIST(1) = NNB
*
*       Perform total force & energy corrections (new polynomial set later).
          KW = KW1
*       Delay velocity kick until routine MDOT on type 13/14/15 in ROCHE.
          IF (IQCOLL.EQ.3.AND.KW1.GE.13) KW = 6
          CALL FCORR(I1,DM,KW)
*
*       Specify commensurate time-step (not needed for BODY(I1) = 0).
          CALL DTCHCK(TIME,STEP(I1),DTK(40))
*
*       Set IPHASE = -3 to preserve ILIST and skip GPSEND in FPOLYI.
          IPHASE = -3
          ISEND = -1
*
*       Initialize new polynomials of neighbours & #I for DM > 0.1 DMSUN.
          IF (DM*ZMBAR.GT.0.1) THEN
*
*       Include body #I at the end (counting from location #2).
              NNB2 = NNB + 2
              ILIST(NNB2) = I1
              NNB2 = NNB2 - 1
*
*       Obtain new F & FDOT and time-steps.
              DO 30 L = 2,NNB2
                  J = ILIST(L)
                  IF (L.EQ.NNB2) THEN
                      J = I1
                  ELSE IF (T0(J).LT.TIME) THEN
                      CALL XVPRED(J,-2)
                      CALL DTCHCK(TIME,STEP(J),DTK(40))
                  END IF
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
                  IF (BODY(J).EQ.0.0D0) THEN
                      CALL FPOLY1(J,J,0)
                  ELSE
                      CALL FPOLYI(J)
                  END IF
   30         CONTINUE
          END IF
          ISEND = -1
          TPREV = TIME - STEPX
      END IF
*
*       See whether closest neighbour forms a KS pair (skip chain).
      IF (IPAIR.GT.0.AND.BODY(I1).GT.0.0D0) THEN
          IF (JMIN.LT.N.AND.RIJ2.LT.RMIN2) THEN
              DO 35 K = 1,3
                  X0DOT(K,JMIN) = XDOT(K,JMIN)
   35         CONTINUE
              ICOMP = MIN(I1,JMIN)
              JCOMP = MAX(I1,JMIN)
              CALL KSREG
              WRITE (6,36) NAME(ICOMP), NAME(JCOMP), LIST(1,2*NPAIRS-1),
     &                     R(NPAIRS), H(NPAIRS), STEP(NTOT)
   36         FORMAT (' COAL KS    NM NP R H DTCM  ',2I6,I4,1P,3E10.2)
              I2 = JMIN
          ELSE
*       Initialize force polynomial for new single body.
              ISEND = -1
              IPHASE = 0
              CALL FPOLYI(ICOMP)
          END IF
      END IF
*
*       Update energy loss & collision counters.
      ECOLL = ECOLL + EB
      E(10) = E(10) + EB
      NPOP(9) = NPOP(9) + 1
      NCOAL = NCOAL + 1
*
*       Open unit #12 the first time.
      IF (FIRST) THEN
          OPEN (UNIT=12,STATUS='NEW',FORM='FORMATTED',FILE='COAL')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF (NCOAL.EQ.1) THEN
              WRITE (12,40)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                       NBIN0, NZERO
   40         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I4,'  N0 =',I6,//)
          END IF
*
          WRITE (12,45)
   45     FORMAT ('   TIME  NAME  NAME  K1  K2  IQ  M1   M2',
     &            '   DM    R1     R2    r/Rc   R     P',/)
      END IF
*
      WRITE (12,50)  TTOT, NAME1, NAME2, KSTAR(I1), KSTAR(I2), IQCOLL,
     &               M1, M2, DM*ZMBAR, RADIUS(I1)*SU, RADIUS(I2)*SU,
     &               RI, RIJ*SU, TK
   50 FORMAT (1X,F6.1,2I6,3I4,3F5.1,2F7.2,F6.2,F7.2,1P,E9.1)
      CALL FLUSH(12)
*
      WRITE (6,55)  WHICH1, IQCOLL, NAME1, NAME2, KSTAR(I1), KSTAR(I2),
     &              KW1, ZMNEW*ZMBAR, RCOLL, EB, DP, DM*ZMBAR, VINF
   55 FORMAT (/,A8,'COAL    IQ =',I3,'  NAME =',2I6,'  K* =',3I3,
     &             '  M =',F5.2,'  RCOLL =',1P,E8.1,'  EB =',E9.1,
     &             '  DP =',E9.1,'  DM =',0P,F6.2,'  V =',F4.1)
*
      KSTAR(I1) = KW1
      KSTAR(I2) = 15
*       Specify IPHASE < 0 for new sorting and ISEND < 0 for new GPSEND.
      IPHASE = -1
      ISEND = -1
      IQCOLL = 0
*
      RETURN
*
      END
