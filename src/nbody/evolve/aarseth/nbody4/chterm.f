      SUBROUTINE CHTERM(ISUB)
*
*
*       Termination of chain system.
*       ----------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,R2(NMX,NMX)
      INTEGER  IJ(NMX)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/ECHAIN/  ECH
*
*
*       Decide between standard termination or collision (ISUB > 0 or < 0).
      IF (ISUB.NE.0) THEN
          ITERM = ISUB
          ISUB = IABS(ISUB)
      END IF
*
*       Prepare KS regularization and direct integration of any other bodies.
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      I5 = 0
      IF (NCH.EQ.2) I3 = I2
      IF (NCH.LE.3) THEN
          I4 = I3
          I5 = I3
          R2(I2,I4) = R2(I2,I3)
          R2(I3,I4) = 0.0
      ELSE IF (NCH.EQ.4) THEN
          I4 = IJ(4)
          I5 = I4
      ELSE IF (NCH.GT.4) THEN
*       Determine indices for second closest pair.
          RX1 = 1.0
          RX0 = R2(I1,I2)
          DO 2 J1 = 1,NCH
*       Avoid choosing close pair I1-I2.
              IF (J1.EQ.I1.OR.J1.EQ.I2) GO TO 2
              DO 1 J2 = J1+1,NCH
              IF (J2.EQ.I1.OR.J2.EQ.I2) GO TO 1
                  IF (R2(J1,J2).LT.RX1.AND.R2(J1,J2).GT.RX0) THEN
                      RX1 = R2(J1,J2)
                      I3 = J1
                      I4 = J2
                  END IF
    1         CONTINUE
    2     CONTINUE
*       Identify remaining single particle(s) by exclusion.
          DO 3 I = 1,NCH
              IF (I.EQ.I1.OR.I.EQ.I2.OR.I.EQ.I3.OR.I.EQ.I4) GO TO 3
              IF (I5.EQ.0) THEN
                  I5 = I
              ELSE
                  I6 = I
              END IF
    3     CONTINUE
      END IF
*
      IF (KZ(30).GT.2) THEN
          WRITE (6,4)  SQRT(R2(I1,I2)), SQRT(R2(I1,I3)),SQRT(R2(I2,I3)),
     &                 SQRT(R2(I2,I4)), SQRT(R2(I3,I4))
    4     FORMAT (' CHTERM:   RIJ (1-2 1-3 2-3 2-4 3-4)  ',1P,5E9.1)
      END IF
*
      JLIST(6) = NAMEC(I1)
      JLIST(7) = NAMEC(I2)
      JLIST(8) = NAMEC(I3)
      JLIST(9) = NAMEC(I4)
      JLIST(10) = NAMEC(I5)
*
*       Specify chain phase indicator and restore original name to c.m. body.
      IPHASE = 8
      NAME(ICH) = NAME0
*
*       Identify current global indices by searching all particles.
      ICM = 0
      DO 10 J = IFIRST,NTOT
          DO 5 L = 1,NCH
              IF (NAME(J).EQ.JLIST(L+5)) THEN
                  JLIST(L) = J
                  IF (BODY(J).GT.0.0D0) ICM = J
              END IF
    5     CONTINUE
   10 CONTINUE
*
*       Modify identification list for special cases NCH = 2 & NCH = 3.
      IF (NCH.EQ.2) THEN
          I3 = I1
          I4 = I2
          JLIST(3) = JLIST(1)
          JLIST(4) = JLIST(2)
      ELSE IF (NCH.EQ.3) THEN
          JLIST(4) = JLIST(3)
      END IF
*
*       Ensure ICOMP < JCOMP for KS regularization.
      ICOMP = MIN(JLIST(1),JLIST(2))
      JCOMP = MAX(JLIST(1),JLIST(2))
*
*       Copy final coordinates & velocities to standard variables.
      LK = 0
      DO 20 L = 1,NCH
          DO 15 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   15     CONTINUE
   20 CONTINUE
*
*       Quantize the elapsed interval within last block-step (factor 1/8).
      TIME2 = T0S(ISUB) + TIMEC - TPREV
      DT8 = (TBLOCK - TPREV)/8.0D0
      IF (DT8.EQ.0.0D0) DT8 = STEP(ICM)/8.0D0
*
*       Adopt the nearest truncated step (at most 8 subdivisions).
      DT2 = TIME2
      IF (TIME2.GT.0.0D0) THEN
          CALL STEPK(DT2,DTN2)
          DTN = NINT(DTN2/DT8)*DT8
      ELSE
*       Choose negative step if pericentre time < TPREV (cf. iteration).
          DT2 = -DT2
          CALL STEPK(DT2,DTN2)
          DTN = -NINT(DTN2/DT8)*DT8
      END IF
*
*       Update time for new polynomial initializations (also for CMBODY).
      TIME = TPREV + DTN
*
*       Set TIMEC for collision to preserve TIME2 after T0S = TIME in SUBSYS.
      IF (ITERM.LT.0.0) TIMEC = T0S(ISUB) + TIMEC - TIME
*
*       Include safety check just in case and save TIME for possible reset.
      TIME = MIN(TBLOCK,TIME)
      TIME0 = TIME
*
*       Predict current X & XDOT for c.m. and neighbours to order F3DOT.
      CALL XVPRED(ICM,-1)
      NNB1 = LISTCM(1) + 1
      DO 25 L = 2,NNB1
*       Note possibility T0(J) = TIME in routine REDUCE would skip on -2.
          J = LISTCM(L)
          CALL XVPRED(J,-1)
   25 CONTINUE
*
*       Copy c.m. coordinates & velocities.
      DO 30 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
   30 CONTINUE
*
*       Set configuration pointers for KS candidates & distant bodies.
      JLIST(6) = I1
      JLIST(7) = I2
      JLIST(8) = I3
      JLIST(9) = I4
      JLIST(10) = I5
*
*       Place new coordinates in the original locations.
      DO 40 L = 1,NCH
          J = JLIST(L)
*       Compare global name & subsystem name to restore the mass & T0.
          DO 32 K = 1,NCH
              IF (NAME(J).EQ.NAMEC(K)) THEN
                  BODY(J) = BODYC(K)
                  T0(J) = TIME
                  CALL DTCHCK(TIME,STEP(J),DTK(40))
              END IF
   32     CONTINUE
*       Transform to global coordinates & velocities using c.m. values.
          LL = JLIST(L+5)
          DO 35 K = 1,3
              X(K,J) = X4(K,LL) + CM(K)
              XDOT(K,J) = XDOT4(K,LL) + CM(K+3)
              X0(K,J) = X(K,J)
              X0DOT(K,J) = XDOT(K,J)
   35     CONTINUE
          ILIST(L+1) = J
   40 CONTINUE
*
*       Set dominant terms in F, FDOT & PHI for NBLIST (#I3 & higher).
      ILIST(1) = NCH
      DO 50 L = 1,NCH
          J = JLIST(L)
          PHI(J) = -VC**2
          CALL FFDOT(J)
   50 CONTINUE
*
*       Update subsystem COMMON variables unless last or only case.
      IF (ISUB.LT.NSUB) THEN
          DO 60 L = ISUB,NSUB
              DO 55 K = 1,6
                  BODYS(K,L) = BODYS(K,L+1)
                  NAMES(K,L) = NAMES(K,L+1)
   55         CONTINUE
              T0S(L) = T0S(L+1)
              TS(L) = TS(L+1)
              STEPS(L) = STEPS(L+1)
              RMAXS(L) = RMAXS(L+1)
              ISYS(L) = ISYS(L+1)
   60     CONTINUE
      END IF
*
*       Check for collision or CE (only needs coordinates & velocities).
      IF (ITERM.LT.0) THEN
          JLIST(1) = ICOMP
          JLIST(2) = JCOMP
*       See whether re-labelling is required (indices I1 - I4 still local).
          IF (R2(I1,I4).LT.R2(I1,I3).OR.R2(I3,I4).LT.R2(I1,I3)) THEN
              IF (R2(I1,I4).LT.R2(I3,I4)) THEN
*       Switch body #I3 & I4 to give new dominant pair I1 & I3.
                  I = JLIST(4)
                  JLIST(4) = JLIST(3)
                  JLIST(3) = I
              ELSE
*       Set JLIST(6) < 0 to denote that body #I3 & I4 will be new KS pair.
                  JLIST(6) = -1
              END IF
          END IF
*       Define collision distance and indicator for routine CMBODY.
          DMINC = RCOLL
          IPHASE = 9
          GO TO 100
      END IF
*
*       Replace ICM in KS perturber lists by all subsystem members.
      CALL NBREST(ICM,NCH,NPAIRS)
*
*       See whether a second binary is present (NCH >= 4 & RB < RMIN).
      KS2 = 0
      IF (NCH.GE.4.AND.R2(I3,I4).LT.RMIN2) THEN
*       Accept second KS pair if well separated (> 2) from smallest binary.
          IF (R2(I3,I4).LT.0.25*MIN(R2(I1,I3),R2(I2,I4))) THEN
              IF (MAX(JLIST(3),JLIST(4)).LE.N) THEN
                  KS2 = 1
              END IF
          END IF
      END IF
*
*       Specify global indices of least dominant bodies (I4 = I3 if NCH = 3).
      I3 = JLIST(3)
      I4 = JLIST(4)
      IF (NCH.EQ.5) I5 = JLIST(5)
*
*       Copy name of two first perturbers in case of exchange in KSREG.
      IF (LISTC(1).GT.0) THEN
          NAMC2 = NAME(LISTC(2))
          IF (LISTC(1).GT.1) THEN
              NAMC3 = NAME(LISTC(3))
          END IF
      END IF
*
*       Save global names of #I3 & I4 for initialization of second KS pair.
      IF (KS2.GT.0) THEN
          NAME3 = NAME(I3)
          NAME4 = NAME(I4)
      END IF
      STEP3 = 0.0
      STEP4 = 0.0
*
*       Initialize force polynomials & time-steps for body #I3 & #I4.
      IF (NCH.GT.2.AND.KS2.EQ.0) THEN
*         CALL FPOLY1(I3,I3,0)
          CALL FPOLYI(I3)
          IF (NCH.GE.4) THEN
*             CALL FPOLY1(I4,I4,0)
              CALL FPOLYI(I4)
              IF (NCH.EQ.5.AND.I5.GT.0) THEN
                  CALL FPOLYI(I5)
              END IF
          END IF
          STEP3 = STEP(I3)
          STEP4 = STEP(I4)
*
*       Check re-initialization of dormant super-hard binary (second pair).
          IF (MAX(I3,I4).GT.N) THEN
              I = MAX(I3,I4)
              CALL RENEW(I)
          END IF
*       Include initialization of 5th body before KS calls change address.
      ELSE IF (NCH.EQ.5.AND.I5.GT.0) THEN
          CALL FPOLYI(I5)
      END IF
*
      IF (NCH.EQ.5.AND.I5.GT.0) THEN
          WRITE (6,70)  I5, NAME(I5), (F(K,I5),K=1,3), STEP(I5)
   70     FORMAT (' CHTERM:    I5 NM F DT  ',2I6,1P,4E10.2)
      END IF
*
*       Perform KS regularization of dominant components (ICOMP < JCOMP).
      IF (JCOMP.LE.N) THEN
          CALL KSREG
*       Restore the time after possible backwards integration in KSPERI.
          TIME = TIME0
      ELSE
*       Initialize components separately and re-activate dormant binary.
*         CALL FPOLY1(JCOMP,JCOMP,0)
          CALL FPOLYI(JCOMP)
          CALL RENEW(JCOMP)
*         CALL FPOLY1(ICOMP,ICOMP,0)
          CALL FPOLYI(ICOMP)
      END IF
*
*       Check for high-velocity ejection with outwards hyperbolic motion.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      IF (NCH.EQ.3) THEN
          DO 74 K = 1,3
              RIJ2 = RIJ2 + (X(K,I3) - X(K,NTOT))**2
              VIJ2 = VIJ2 + (XDOT(K,I3) - XDOT(K,NTOT))**2
              RDOT = RDOT + (X(K,I3) - X(K,NTOT))*
     &                      (XDOT(K,I3) - XDOT(K,NTOT))
   74     CONTINUE
          HI = 0.5*VIJ2 - (BODY(I3) + BODY(NTOT))/SQRT(RIJ2)
          IF (HI.GT.0.0.AND.RDOT.GT.0.0) THEN
               CALL HIVEL(I3)
          END IF
      END IF
*
      IF (KZ(30).GT.1) THEN
          WRITE (6,75) TIME+TOFF, I3, I4, NTOT, STEP3, STEP4, STEP(NTOT)
   75     FORMAT (' CHTERM:   T I3 I4 NT DT ',F9.4,3I6,1P,3E9.1)
      END IF
*
*       Initialize second KS pair if separation is small.
      IF (KS2.GT.0) THEN
          I3 = 0
*       Re-determine the global indices after updating in KSREG.
          DO 80 I = IFIRST,N
              IF (NAME(I).EQ.NAME3) I3 = I
              IF (NAME(I).EQ.NAME4) THEN
                  I4 = I
                  IF (I3.GT.0) GO TO 85
              END IF
   80     CONTINUE
*
*       Define components and perform new regularization.
   85     ICOMP = MIN(I3,I4)
          JCOMP = MAX(I3,I4)
          CALL KSREG
          TIME = TIME0
*
          IF (KZ(30).GT.1) THEN
              WRITE (6,90)  I3, I4, STEP(NTOT), STEP(2*NPAIRS-1),
     &                      R(NPAIRS), H(NPAIRS), GAMMA(NPAIRS)
   90         FORMAT (' CHTERM:   SECOND BINARY   I3 I4 DT DTK R H G ',
     &                                            2I5,1P,5E9.1)
          END IF
      END IF
*
*       Rename any old perturbers that have been exchanged in KSREG.
      IF (LISTC(1).GT.0.AND.LISTC(2).LT.IFIRST) THEN
          DO 92 J = IFIRST,N
              IF (NAME(J).EQ.NAMC2) THEN
                  LISTC(2) = J
              ELSE IF (LISTC(1).GT.1.AND.NAME(J).EQ.NAMC3) THEN
                  LISTC(3) = J
              END IF
   92     CONTINUE
      END IF
*
*       Initialize the perturbers to improve derivatives (same locations).
      NBC1 = LISTC(1) + 1
      DO 95 L = 2,NBC1
          J = LISTC(L)
          IF (J.GT.N) GO TO 95
          DO 94 K = 1,3
              X0DOT(K,J) = XDOT(K,J)
   94     CONTINUE
          CALL FPOLYI(J)
   95 CONTINUE
*
*       Check minimum two-body distance.
      DMINC = MIN(DMINC,RCOLL)
*
*       Update net binary energy change.
      CHCOLL = CHCOLL + CM(9)
*
*       Update number of DIFSY calls, tidal dissipations & collision energy.
      NSTEPC = NSTEPC + NSTEP1
      NDISS = NDISS + NDISS1
      NSYNC = NSYNC + ISYNC
      ECOLL = ECOLL + ECOLL1
      E(10) = E(10) + ECOLL1
*
*       Reduce subsystem counter and initialize membership & internal energy.
      NSUB = NSUB - 1
      NCH = 0
      ECH = 0.0
*
*       Check for subsystem at last COMMON dump (no restart with NSUB > 0).
      IF (NSUB.EQ.0.AND.KZ(2).GT.1) THEN
          IF (TIME - TDUMP.LT.TIMEC) THEN
              TDUMP = TIME
              CALL MYDUMP(1,2)
          END IF
      END IF
*
*       Specify non-zero indicator for sending new data to GRAPE in INTGRT.
  100 ISEND = -1
*       Set new value of TPREV (note: must be commensurate for any 2nd call).
      IF (ITERM.GE.0) TPREV = TIME - 16.0*DT8
*
      RETURN
*
      END
