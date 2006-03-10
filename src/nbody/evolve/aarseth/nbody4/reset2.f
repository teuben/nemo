      SUBROUTINE RESET2
*
*
*      Termination of double hierarchy. 
*      --------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  XX(3,3),VV(3,3)
      CHARACTER*11  WHICH1
      INTEGER JSAVE(LMAX)
*
*
*       Set index of disrupted pair and save output parameters.
      IPAIR = KSPAIR
      I = N + IPAIR
      E1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
      G1 = GAMMA(IPAIR)
      R1 = R(IPAIR)
      SEMI1 = -0.5*BODY(I)/H(IPAIR)
*
*       Locate current position in the merger table.
      IMERGE = 0
      DO 1 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(I)) IMERGE = K
    1 CONTINUE
*
*       Produce info on quintuplet or sextuplet (NAMEG has outer c.m. name).
      IF (NAMEG(IMERGE).GT.NZERO) THEN
          RI2 = 0.0
          DO 2 K = 1,3
              RI2 = RI2 + (X(K,I) - RDENS(K))**2
    2     CONTINUE
          RI = SQRT(RI2)
          ECC2 = (1.0 - R1/SEMI1)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI1)
          ECC = SQRT(ECC2)
          PMIN = SEMI1*(1.0 - ECC)
          EB = 0.5*BODY(2*IPAIR-1)*BODY(2*IPAIR)/SEMI1
          EB = EB*FLOAT(N-NPAIRS)/ZKIN
          ZM = BODY(I)*SMU
          WHICH1 = ' QUINTUPLET'
*       Determine previous merger index with NAME greater by 2*NZERO.
          JM = 1
          DO 3 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I) + 2*NZERO) JM = K
    3     CONTINUE
*       Use previous merger index to identify ghost from earlier level.
          DO 4 J = IFIRST,NTOT
              IF (NAME(J).EQ.NAMEG(JM)) JG = J
    4     CONTINUE
          IF (JG.GT.N) WHICH1 = ' SEXTUPLET '
          WRITE (6,5)  WHICH1, TTOT, ZM, NAME(2*IPAIR-1), NAME(2*IPAIR), 
     &                 NAMEG(IMERGE), RI, ECC, EB, SEMI1, PMIN, G1
    5     FORMAT (/,' END',A11,'   T MT NM1 NM2 NM3 RI E1 EB1 A1 PM G1',
     &                             F10.2,F6.2,3I6,2F6.2,F6.1,1P,3E10.2)
      END IF
*
*       Check presence of [[B,B],S] quintuplet.
      IF (NAMEG(IMERGE).GT.0.AND.NAMEG(IMERGE).LE.NZERO) THEN
          I1 = 2*IPAIR - 1
          I2 = 2*IPAIR
          CALL FINDJ(I1,JG,IM)
          IF (NAME(JG).GT.NZERO) THEN
              ZM = BODY(I)*SMU
              WRITE (6,9)  TTOT, ZM, NAME(2*IPAIR), NAME(JG), NAME(I2)
    9         FORMAT (/,' END QUINT2    T NM1 NMG NM3 ',F9.2,F6.2,3I6)
          END IF
      END IF
*
*       Include diagnostics for double triple.
      IF (NAMEG(IMERGE).LT.0) THEN
          I1 = 2*IPAIR - 1
          J1 = 2*IPAIR
          CALL FINDJ(I1,JI,IM)
          CALL FINDJ(J1,JJ,JM)
          AI = -0.5*(CM(1,IM) + CM(2,IM))/HM(IM)
          AJ = -0.5*(CM(1,JM) + CM(2,JM))/HM(JM)
          ECC2 = (1.0 - R1/SEMI1)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI1)
          ECC = SQRT(ECC2)
          PMIN = SEMI1*(1.0 - ECC)
          ZM = BODY(I)*SMU
          WRITE (6,11)  TTOT, NAMEG(IM), NAME(JI), NAMEG(JM), NAME(JJ),
     &                  AI, AJ, R(IPAIR), SEMI1, PMIN
   11     FORMAT (/,' END HITRIP    T MT NM AI AJ R1 A1 PM ',
     &                              F9.2,4I6,1P,5E10.2)
      END IF
*
*       Check optional diagnostics.
      IF (KZ(18).GT.0.AND.KSTARM(IMERGE).LE.10) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Save perturbers for correction procedure.
      I1 = 2*IPAIR - 1
      NNB = LIST(1,I1) + 1
      DO 6 L = 2,NNB
          J = LIST(L,I1)
          JPERT(L) = J
    6 CONTINUE
*
*       Ensure that c.m. coordinates are known to highest order.
      CALL XVPRED(I,0)
*
*       Predict perturber coordinates & velocities (XDOT used by FPOLY1).
      DO 7 L = 2,NNB
          J = JPERT(L)
          CALL XVPRED(J,0)
    7 CONTINUE
*
*       Obtain current coordinates & velocities and specify KS components.
      CALL RESOLV(IPAIR,2)
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
*
*       Initialize mass, coordinates and velocities for new cm body.
      BODY(I) = BODY(ICOMP)
      DO 8 K = 1,3
          X(K,I) = X(K,ICOMP)
          XDOT(K,I) = XDOT(K,ICOMP)
          X0DOT(K,I) = XDOT(K,ICOMP)
    8 CONTINUE
*
*       Add outer component to perturber list.
      JPERT(1) = JCOMP
*
*       Sum first part of potential energy correction due to tidal effect.
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT1)
*
*       Find correct location of ghost particle using identification name.
      ICM = I
      JCOMP1 = JCOMP
      DO 10 I = 1,NTOT
          IF (BODY(I).EQ.0.0D0.AND.NAME(I).EQ.NAMEG(IMERGE)) JCOMP1 = I
   10 CONTINUE
*
*       Regularize two-body configuration if JCOMP1 cannot be identified.
      IF (JCOMP.EQ.JCOMP1) THEN
          WRITE (6,12)  IMERGE, NAMEG(IMERGE), JCOMP
   12     FORMAT (/,5X,'WARNING!    RESET2    JCOMP NOT IDENTIFIED ',
     &                       '   IM =',I3,'  NAMEG =',I6,'  JCOMP =',I6)
          GO TO 100
      END IF
*
*       Initialize basic variables for ghost and new c.m (JCOMP -> JCOMP1).
      J1 = JCOMP1
      J = JCOMP
   13 T0(J1) = TIME
      BODY(J1) = BODY(J)
      DO 14 K = 1,3
          X(K,J1) = X(K,J)
          X0(K,J1) = X(K,J)
          XDOT(K,J1) = XDOT(K,J)
          X0DOT(K,J1) = XDOT(K,J)
   14 CONTINUE
      IF (J.EQ.JCOMP) THEN
          J1 = ICM
          J = ICOMP
          GO TO 13
      END IF
*
*       Restore masses, coordinates & velocities of old hierarchical binary.
      BODY(ICOMP) = CM(1,IMERGE)
      BODY(JCOMP) = CM(2,IMERGE)
      ZM = -BODY(ICOMP)/(BODY(ICOMP) + BODY(JCOMP))
*
*       Begin with second component since ICOMP holds new c.m. variables.
      I = JCOMP
      DO 20 KCOMP = 1,2
          DO 15 K = 1,3
              X(K,I) = X(K,ICOMP) + ZM*XREL(K,IMERGE)
              X0DOT(K,I) = X0DOT(K,ICOMP) + ZM*VREL(K,IMERGE)
              XDOT(K,I) = X0DOT(K,I)
*       Note that XDOT is only needed for improved polynomials of JCOMP.
   15     CONTINUE
          I = ICOMP
          ZM = BODY(JCOMP)/(BODY(ICOMP) + BODY(JCOMP))
   20 CONTINUE
*
*       Copy KS variables for inner binary (small TDOT2 near apo/peri).
      T0(I1) = TIME
      LIST(1,I1) = 1
      H(IPAIR) = HM(IMERGE)
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
      VI2 = 0.0
      DO 30 K = 1,4
          U(K,IPAIR) = UM(K,IMERGE)
          U0(K,IPAIR) = U(K,IPAIR)
          UDOT(K,IPAIR) = UMDOT(K,IMERGE)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          VI2 = VI2 + UDOT(K,IPAIR)**2
   30 CONTINUE
*
*       Include safety check to identify possible bug of January 95.
      HI = (2.0*VI2 - (BODY(ICOMP) + BODY(JCOMP)))/R(IPAIR)
      IF (ABS(HI - H(IPAIR)).GT.0.01*ABS(HI)) THEN
          TM = MIN(TEV(ICOMP),TEV(JCOMP))
          WRITE (6,35)  TTOT, TM, NMERGE, IMERGE, R(IPAIR), H(IPAIR), HI
   35     FORMAT (' DANGER!    RESET2    T TM NM IM R H0 HI ',
     &                                   2F9.2,2I4,1P,3E10.2)
          H(IPAIR) = HI
      END IF
*
*       Initialize force polynomial for outer component using resolved c.m.
      CALL FPOLY1(JCOMP1,JCOMP1,0)
*
*       See whether body #JCOMP1 should be included in NLIST.
      IF (T0(JCOMP1) + STEP(JCOMP1).LT.TLIST) THEN
          CALL NLMOD(JCOMP1,1)
      END IF
*
*       Save perturber list to prevent over-writing in FPOLYI (FPOLY1 is OK).
      DO 40 L = 1,NNB
          JSAVE(L) = JPERT(L)
   40 CONTINUE
*
*       Initialize c.m. polynomials and activate inner binary.
      CALL KSIN2(3)
*
*       Restore original name of inner hierarchy (c.m. NAME set in MERGE2).
      NAME(ICM) = NAME(ICM) + 2*NZERO
*
*       Locate current position in the merger table (IM = 1 for safety).
      IM = 1
      DO 42 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(ICM)) IM = K
   42 CONTINUE
*
*       Evaluate stability parameter from current elements.
      Q = BODY(JCOMP)/BODY(ICOMP)
      SEMI0 = -0.5*BODY(ICOMP)/HM(IM)
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      ECC1 = SQRT(ECC2)
      XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      RSTAB = 2.8*XFAC**0.4*SEMI0
*
*       Determine inclination between inner relative motion and outer orbit.
      DO 43 K = 1,3
          XX(K,1) = XREL(K,IM)
          XX(K,2) = 0.0
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = VREL(K,IM)
          VV(K,2) = 0.0
          VV(K,3) = XDOT(K,JCOMP)
   43 CONTINUE
      CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ALPHA)
*
*       Include standard inclination factor in the stability criterion.
      YFAC = 1.0 - 0.3*ALPHA/180.0
      R0(IPAIR) = YFAC*RSTAB
*
*       Restore perturber list for correction procedure.
      DO 44 L = 1,NNB
          JPERT(L) = JSAVE(L)
   44 CONTINUE
*
*       Rename perturber list for routine NBPOT.
      JPERT(1) = JCOMP1
*
*       Restore Roche stage indicator and binary flag (any ghost c.m. is OK).
      KSTAR(ICM) = KSTARM(IMERGE)
      LIST(2,JCOMP) = IFLAGM(IMERGE)
*
*       See whether the outer component is a single or composite particle.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Restore component masses for outer binary.
      JPAIR = JCOMP1 - N
      BODY(2*JPAIR-1) = CM(3,IMERGE)
      BODY(2*JPAIR) = CM(4,IMERGE)
*
*       Update look-up times & radii and check possible Roche condition.
      IF (KZ(19).GE.3.AND.KZ(34).GT.0) THEN
          IF (KSTAR(JCOMP1).GT.0.AND.KSTAR(JCOMP1).LE.20) THEN
              CALL TRFLOW(JPAIR,DTR)
              TEV(JCOMP1) = TIME + DTR
              TMDOT = MIN(TEV(JCOMP1),TMDOT)
              TMDOT = MIN(TEV(2*JPAIR),TMDOT)
          END IF
      END IF
*
*       Obtain coordinates & velocities of unperturbed binary components.
      CALL RESOLV(JPAIR,1)
*
*       Select new perturbers and initialize polynomials for KS motion.
      CALL KSLIST(JPAIR)
      CALL KSPOLY(JPAIR,1)
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
      CALL NBPOT(2,NNB,POT3)
      JLIST(1) = JCOMP1
      CALL NBPOT(1,NNB,POT4)
*
*       Update the merger energy.
      EB2 = BODY(2*JPAIR-1)*BODY(2*JPAIR)*H(JPAIR)/BODY(JCOMP1)
      EMERGE = EMERGE - EB2
*
      E2 = E1/EB2
      EB2 = EB2/BE(3)
      DP = POT4 - POT3
      IF (KZ(15).GT.1.OR.G1.GT.GMAX) THEN
          WRITE (6,45)  JPAIR, H(JPAIR), BODY(2*JPAIR-1),
     &                  BODY(2*JPAIR), R1, SEMI1, G1, E2, EB2,
     &                  R(JPAIR), GAMMA(JPAIR), DP
   45     FORMAT (' END OUTER QUAD',I4,'  H =',F7.1,'  M =',2F7.4,

     &            '  R1 =',1P,E8.1,'  A1 =',E8.1,'  G1 =',0P,E8.1,
     &            '  E1 =',F6.3,'  EB2 =',F6.3,'  RB2 =',1P,E8.1,
     &            '  G2 =',E8.1,'  DP =',E8.1)
      END IF
*
*       Include interaction of body #ICOMP & JCOMP with perturbers.
   50 JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT2)
*
*       Form square of c.m. velocity correction due to tidal effects.
      VI2 = X0DOT(1,ICM)**2 + X0DOT(2,ICM)**2 + X0DOT(3,ICM)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
      CORR = 1.0 + 2.0*DPHI/(BODY(ICM)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
*     DO 60 K = 1,3
*         X0DOT(K,ICM) = SQRT(CORR)*X0DOT(K,ICM)
*  60 CONTINUE
*
*       Modify the merger energy to maintain conservation.
      EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(ICM)
      EMERGE = EMERGE - EB + DPHI
*
      E1 = E1/EB
      EB = EB/BE(3)
      IF (KZ(15).GT.1) THEN
          WRITE (6,65)  IMERGE, TTOT, BODY(2*IPAIR-1),
     &                  BODY(2*IPAIR), R1, SEMI1, EB, E1,
     &                  GAMMA(IPAIR), G1, NNB-1
   65     FORMAT (' END MERGE2',I3,'  T =',F8.2,'  M =',2F7.4,
     &            '  R1 =',1PE8.1,'  A1 =',E8.1,'  EB =',0PF6.3,
     &            '  E1 =',F6.3,'  GB =',1PE8.1,'  G =',0PF6.3,
     &            '  NB =',I3)
      END IF
*
*       Check Roche look-up time.
      IF (KSTAR(ICM).GT.0.AND.KSTAR(ICM).LE.20.AND.KZ(34).GT.0) THEN
          K = ICM - N
          CALL TRFLOW(K,DTR)
          TEV(ICM) = MIN(TEV(ICM),TIME + DTR)
          TMDOT = MIN(TEV(ICM),TMDOT)
      END IF
*
*       Reduce merger counter and update tables (unless last or only pair).
   70 NMERGE = NMERGE - 1
      DO 80 L = IMERGE,NMERGE
          L1 = L + 1
          HM(L) = HM(L1)
          TMDIS(L) = TMDIS(L1)
          NAMEG(L) = NAMEG(L1)
          NAMEM(L) = NAMEM(L1)
          KSTARM(L) = KSTARM(L1)
          IFLAGM(L) = IFLAGM(L1)
          DO 74 K = 1,3
              XREL(K,L) = XREL(K,L1)
              VREL(K,L) = VREL(K,L1)
   74     CONTINUE
          DO 75 K = 1,4
              CM(K,L) = CM(K,L1)
              UM(K,L) = UM(K,L1)
              UMDOT(K,L) = UMDOT(K,L1)
   75     CONTINUE
   80 CONTINUE
*
*       Examine merger list for possible escapers (retain up to 3 levels).
      DO 90 L = 1,NMERGE
          DO 85 J = 1,NPAIRS
              IF (NAMEM(L).EQ.NAME(N+J).OR.
     &            NAMEM(L).EQ.NAME(N+J) + 2*NZERO.OR.
     &            NAMEM(L).EQ.NAME(N+J) + 4*NZERO) GO TO 90
   85     CONTINUE
*       Remove tables for any merger not identified.
          IMERGE = L
          GO TO 70
   90 CONTINUE
*
*       Specify non-zero indicator for sending new data to GRAPE in INTGRT.
      ISEND = -1
*       Set IPHASE < 0 for new sorting.
      IPHASE = -1
*
  100 RETURN
*
      END
