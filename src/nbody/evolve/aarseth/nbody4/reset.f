      SUBROUTINE RESET
*
*
*       Restore hierarchical configuration.
*       -----------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      INTEGER JSAVE(5*LMAX)
*
*
*       Set indices of disrupted pair.
      IPAIR = KSPAIR
      I = N + IPAIR
*
*       Check special termination of double hierarchy (routine RESET2).
      IF (NAME(I).LT.-2*NZERO) THEN
          CALL RESET2
          GO TO 100
      END IF
*
*       Check termination of double hierarchy identified by zero mass.
      IF (BODY(I).EQ.0.0D0) THEN
          DO 2 K = 1,NMERGE
              IF (NAME(I).EQ.NAMEG(K)) THEN
                  KK = 0
*       Note possible termination of the wrong hierarchy if more than one.
                  DO 1 J = N+1,NTOT
                      IF (NAME(J).LT.-2*NZERO) KK = J
    1             CONTINUE
                  IF (KK.EQ.0) THEN
                      WRITE (6,9)  I, NAME(I)
    9                 FORMAT (' DANGER!    RESET    I NM  ',2I6)
                      CALL gpfree
                      STOP
                  END IF
                  IF (KK.GT.0) KSPAIR = KK - N
                  CALL RESET2
                  GO TO 100
              END IF
    2     CONTINUE
      END IF
*
*       Save diagnostic parameters.
      E1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
      G1 = GAMMA(IPAIR)
      R1 = R(IPAIR)
      SEMI1 = -0.5*BODY(I)/H(IPAIR)
*
*       Locate current position in the merger table.
      IMERGE = 0
      DO 3 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(I)) IMERGE = K
    3 CONTINUE
*
*       Check optional diagnostics.
      IF (KZ(18).GT.0.AND.KSTARM(IMERGE).LE.10) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Save perturbers for correction procedure and rename if moved up.
      I1 = 2*IPAIR - 1
      NNB = LIST(1,I1) + 1
      DO 4 L = 2,NNB
          J = LIST(L,I1)
          IF (J.GT.I) J = J - 1
          IF (J.LE.2*NPAIRS.AND.J.GT.2*IPAIR) J = J - 2
          JSAVE(L) = J
    4 CONTINUE
*
*       Ensure that c.m. coordinates are known to highest order.
      CALL XVPRED(I,0)
*
*       Restore original c.m. name and terminate outer KS pair.
      NAME(I) = NZERO - NAME(I)
      CALL KSTERM
*
*       Predict perturber coordinates & velocities (XDOT used by FPOLY1).
      DO 5 L = 2,NNB
          JPERT(L) = JSAVE(L)
          J = JPERT(L)
          CALL XVPRED(J,0)
    5 CONTINUE
*
*       Add outer component to perturber list and set old dipole range.
      JPERT(1) = JCOMP
*     RCRIT2 = CMSEP2*(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
*    &                                                XREL(3,IMERGE)**2)
*
*       Sum first part of potential energy correction due to tidal effect.
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT1)
*
*       Find the nearest neighbour and reduce steps of active perturbers.
*     RJMIN2 = 100.0
*     JMIN = N
*     DO 8 L = 1,NNB
*         J = JPERT(L)
*         RIJ2 = (X(1,J) - X(1,ICOMP))**2 + (X(2,J) - X(2,ICOMP))**2 +
*    &                                      (X(3,J) - X(3,ICOMP))**2
*       Identify the closest perturber for dominant force calculation.
*         IF (RIJ2.LT.RJMIN2.AND.J.NE.JCOMP) THEN
*             RJMIN2 = RIJ2
*             JMIN = J
*         END IF
*         IF (RIJ2.LT.RCRIT2) THEN
*       Reduce step of inner binary perturbers (c.m. approximation used).
*             STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*       Include particle #J in time-step list unless already present.
*             IF (T0(J) + STEP(J).LT.TLIST) THEN
*                 CALL NLMOD(J,1)
*             END IF
*         END IF
    8 CONTINUE
*
*       Find correct location of ghost particle using identification name.
      JCOMP1 = JCOMP
*       Note that ghost may be saved in an old binary c.m. (search NTOT). 
      DO 10 I = 1,NTOT
          IF (BODY(I).EQ.0.0D0.AND.NAME(I).EQ.NAMEG(IMERGE)) JCOMP = I
   10 CONTINUE
*
*       Regularize two-body configuration if JCOMP cannot be identified.
      IF (JCOMP.EQ.JCOMP1) THEN
          WRITE (6,12)  NAMEG(IMERGE)
   12     FORMAT (/,5X,'DANGER!    RESET    JCOMP NOT IDENTIFIED ',
     &                                                  '   NAMEG =',I5)
          CALL gpfree
          STOP
*         CALL KSREG
*         GO TO 70
      END IF
*
*       Restore masses, coordinates & velocities of inner binary.
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
*       Add JCOMP to perturber lists containing ICOMP (KSREG sets c.m.).
      JLIST(1) = JCOMP
      CALL NBREST(ICOMP,1,NPAIRS)
*
*       Initialize force polynomial for outer component using resolved c.m.
      CALL FPOLY1(JCOMP1,JCOMP1,0)
*
*       See whether body #JCOMP1 should be included in NLIST.
      IF (T0(JCOMP1) + STEP(JCOMP1).LT.TLIST) THEN
          CALL NLMOD(JCOMP1,1)
      END IF
*
*       Rectify merger variables before copying.
      CALL HIRECT(IMERGE)
*
*       Copy basic KS variables for inner binary (small TDOT2 near apo/peri).
      JP1 = NPAIRS + 1
      H(JP1) = HM(IMERGE)
      R(JP1) = 0.0D0
      TDOT2(JP1) = 0.0D0
      VI2 = 0.0
      DO 30 K = 1,4
          U(K,JP1) = UM(K,IMERGE)
          U0(K,JP1) = U(K,JP1)
          UDOT(K,JP1) = UMDOT(K,IMERGE)
          R(JP1) = R(JP1) + U(K,JP1)**2
          VI2 = VI2 + UDOT(K,JP1)**2
   30 CONTINUE
*
*       Include safety check to identify possible bug of January 95.
      HI = (2.0*VI2 - (BODY(ICOMP) + BODY(JCOMP)))/R(JP1)
      IF (ABS(HI - H(JP1)).GT.0.01*ABS(HI)) THEN
          TM = MIN(TEV(ICOMP),TEV(JCOMP))
          WRITE (6,35)  TTOT, TM+TOFF, NMERGE, IMERGE, R(JP1), H(JP1), HI
   35     FORMAT (' DANGER!    RESET    T TM NM IM R H0 HI ',
     &                                  2F9.2,2I4,1P,3E10.2)
*         H(JP1) = HI
      END IF
*
*       Rename perturber list for routine NBPOT.
      JPERT(1) = JCOMP
*
*       Save perturber list to prevent over-writing in FPOLYI (FPOLY1 is OK).
      DO 40 L = 1,NNB
          JSAVE(L) = JPERT(L)
   40 CONTINUE
*
*       Save ghost index and re-activate inner binary (JCOMP <-> JCOMP1).
      JCOMP1 = JCOMP
      CALL KSREG
*
*       Restore perturber list for correction procedure.
      DO 44 L = 1,NNB
          JPERT(L) = JSAVE(L)
   44 CONTINUE
*
*       Restore Roche stage indicator and binary flag (any ghost c.m. is OK).
      KSTAR(N+NPAIRS) = KSTARM(IMERGE)
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
              TMDOT = MIN(TEV(JCOMP1),TEV(2*JPAIR),TMDOT)
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
      IF (KZ(15).GT.1) THEN
          WRITE (6,45)  JPAIR, H(JPAIR), BODY(2*JPAIR-1),
     &                  BODY(2*JPAIR), R1, SEMI1, G1, E2, EB2,
     &                  R(JPAIR), GAMMA(JPAIR), DP
   45     FORMAT (' END QUAD',I4,'  H =',F7.1,'  M =',2F7.4,
     &            '  R1 = ',1P,E8.1,'  A1 =',E8.1,'  G1 =',0P,F7.3,
     &            '  E2 =',F6.3,'  EB2 =',F6.3,'  RB2 =',1PE8.1,
     &            '  G2 =',E8.1,'  DP =',E8.1)
      END IF
*
*       Include interaction of body #ICOMP & JCOMP with perturbers.
   50 JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT2)
*
*       Form square of c.m. velocity correction due to tidal effects.
      VI2 = X0DOT(1,NTOT)**2 + X0DOT(2,NTOT)**2 + X0DOT(3,NTOT)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
*     CORR = 1.0 + 2.0*DPHI/(BODY(NTOT)*VI2)
*     IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
*     DO 60 K = 1,3
*         X0DOT(K,NTOT) = SQRT(CORR)*X0DOT(K,NTOT)
*  60 CONTINUE
*
*       Modify the merger energy to maintain conservation.
      EB = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)*H(NPAIRS)/BODY(NTOT)
*       Note that EMERGE may contain escaped mergers.
      EMERGE = EMERGE - EB + DPHI
*
      E1 = E1/EB
      EB = EB/BE(3)
      IF (KZ(15).GT.1) THEN
          WRITE (6,65)  IMERGE, TTOT, BODY(2*NPAIRS-1),
     &                  BODY(2*NPAIRS), R1, SEMI1, EB, E1,
     &                  GAMMA(NPAIRS), G1, NNB-1
   65     FORMAT (' END MERGER',I3,'  T =',F8.2,'  M =',2F7.4,
     &            '  R1 =',1PE8.1,'  A1 =',E8.1,'  EB =',0PF6.3,
     &            '  E1 =',F6.3,'  GB =',1PE8.1,'  G =',0PF6.3,
     &            '  NB =',I3)
      END IF
*
*       Check Roche look-up time for circular orbit and ensure priority.
      IF (KSTAR(NTOT).GE.10.AND.KSTAR(NTOT).LE.20) THEN
          CALL TRFLOW(NPAIRS,DTR)
          TEV(NTOT) = TIME + DTR
          TMDOT = MIN(TEV(NTOT),TEV(ICOMP),TMDOT)
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
*       Include consistency check on merger names.
      DO 96 I = N+1,NTOT
          IF (NAME(I).LT.0) THEN
              DO 94 L = 1,NMERGE
                  IF (NAMEM(L).EQ.NAME(I)) GO TO 96
   94         CONTINUE
              WRITE (6,95)  I, NAME(I), (NAMEM(L),L=1,NMERGE)
   95         FORMAT (' DANGER!    RESET    I NAMI NAMEM  ',10I6)
              CALL gpfree
              STOP
          END IF
   96 CONTINUE
*
  100 RETURN
*
      END
