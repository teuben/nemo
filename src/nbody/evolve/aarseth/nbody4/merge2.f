      SUBROUTINE MERGE2
*
*
*       Merging of double hierarchy.
*       ----------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      CHARACTER*11  WHICH1
*
*
*       COMMON variables for merge procedure
*       ************************************
*
*       ------------------------------------------------------------------
*       CM      Component masses of merged binary (2nd binary in CM(3&4)).
*       HM      Binding energy per unit mass of merged binary.
*       IFLAGM  Flag of merged binary (<0: primordial; =0: new binary).
*       KSTARM  Roche stage indicator (KSTAR for 1st c.m.; ghost is OK).
*       NAMEG   Name of the associated ghost component.
*       NAMEM   Name of new c.m. (< 0) for identification of merger index.
*       NMERG   Total number of mergers.
*       NMERGE  Current number of merged binaries (maximum is MMAX).
*       UM      Regularized coordinates of merged binary.
*       UMDOT   Regularized velocity of merged binary.
*       VREL    Relative velocity of merged binary components.
*       XREL    Relative coordinates of merged binary components.
*       ------------------------------------------------------------------
*
*
*       Ensure that the hierarchy is retained as primary KS pair.
      IF (NAME(N+KSPAIR).GT.0) THEN
          K = KSPAIR
          KSPAIR = JCOMP - N
          JCOMP = N + K
      END IF
*
*       Set pair index & c.m. of inner binary and save outer component.
      IPAIR = KSPAIR
      I = N + IPAIR
      JCOMP1 = JCOMP
      ICOMP1 = I
*
*       Produce diagnostics for standard quintuplet or sextuplet system.
      IF (NAME(JCOMP).GT.NZERO) THEN
          SEMI = -0.5*BODY(I)/H(IPAIR)
          RI2 = 0.0
          DO 2 K = 1,3
              RI2 = RI2 + (X(K,I) - RDENS(K))**2
    2     CONTINUE
          E2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(E2)
          RI = SQRT(RI2)
          EB = 0.5*BODY(2*IPAIR-1)*BODY(2*IPAIR)/SEMI
          EB = EB*FLOAT(N-NPAIRS)/ZKIN
          EB = MIN(EB,999.9D0)
          ZM = (BODY(I) + BODY(JCOMP))*SMU
          WHICH1 = ' QUINTUPLET'
          CALL FINDJ(I,JG,IM)
          IF (JG.GT.N) WHICH1 = ' SEXTUPLET '
          WRITE (6,5)  WHICH1, TTOT, ZM, NAME(2*IPAIR-1), NAME(2*IPAIR),
     &                 NAME(JCOMP), RI, ECC, EB, SEMI,PCRIT,GAMMA(IPAIR)
    5     FORMAT (/,' NEW',A11,'   T MT NM1 NM2 NM3 RI E0 EB0 A0 PC G0',
     &                             F10.2,F6.2,3I6,2F6.2,F6.1,1P,3E10.2)
      END IF
*
*       Check formation of [[B,B],S] quintuplet (expected to be rare).
      IF (NAME(JCOMP).GT.0.AND.NAME(JCOMP).LE.NZERO) THEN
          CALL FINDJ(I,JG,IM)
          IF (NAME(JG).GT.NZERO) THEN
              ZM = (BODY(I) + BODY(JCOMP))*SMU
              WRITE (6,8)  TTOT, ZM, NAME(2*IPAIR), NAME(JG),NAME(JCOMP)
    8         FORMAT (/,' NEW QUINT2    T MT NM1 NMG NM3 ',
     &                                  F9.2,F6.2,3I6)
          END IF
      END IF
*
*       Include diagnostics for double triple ([[B,S],[B,S]]).
      IF (NAME(JCOMP).LT.0) THEN
          JPAIR = JCOMP - N
          I1 = 2*IPAIR - 1
          J1 = 2*JPAIR - 1
          CALL FINDJ(I1,JI,IM)
          CALL FINDJ(J1,JJ,JM)
          AI = -0.5*(CM(1,IM) + CM(2,IM))/HM(IM)
          AJ = -0.5*(CM(1,JM) + CM(2,JM))/HM(JM)
          ZM = (BODY(I) + BODY(JCOMP))*SMU
          GX = MAX(GAMMA(IPAIR),GAMMA(JPAIR))
          WRITE (6,10)  TTOT, ZM, NAME(I1), NAME(2*IPAIR), NAME(J1),
     &                  NAME(2*JPAIR), AI, AJ, R(IPAIR), R(JPAIR),
     &                  PCRIT, GX
   10     FORMAT (/,' NEW HITRIP    T MT NM AI AJ RI RJ PC GX ',
     &                              F9.2,F6.2,4I6,1P,6E10.2)
      END IF
*
*       Ensure correct coordinates & velocities.
      CALL RESOLV(IPAIR,3)
*
*       Check optional diagnostics.
      IF (KZ(18).GT.0.AND.KSTAR(I).LE.10) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Increase merger counter and set current merger index.
      NMERG = NMERG + 1
      NMERGE = NMERGE + 1
      IMERGE = NMERGE
*
*       Save component masses and evaluate reduced mass of inner binary.
      CM(1,IMERGE) = BODY(2*IPAIR-1)
      CM(2,IMERGE) = BODY(2*IPAIR)
      CM(3,IMERGE) = 0.0D0
      ZMU = BODY(2*IPAIR-1)*BODY(2*IPAIR)/BODY(N+IPAIR)
*
*       Set current energy of any second perturbed binary (from XVPRED).
      IF (JCOMP.GT.N) THEN 
          CALL XVPRED(JCOMP,-1)
          JPAIR = JCOMP - N
          IF (LIST(1,2*JPAIR-1).GT.0) H(JPAIR) = HT
*       Ensure that outer binary components are resolved.
          CALL KSRES(JPAIR,J1,J2,0.0D0)
*       Enforce non-zero membership for potential correction in NBPOT.
          LIST(1,2*JPAIR-1) = 1
      ELSE
          CALL XVPRED(JCOMP,-1)
      END IF
*
*       Update the primary velocity of body #JCOMP.
      DO 15 K = 1,3
          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
   15 CONTINUE
*
*       Include more perturbers in tidal correction (maybe not important).
*     RS2 = CMSEP2*RMIN2
*     IP = IPHASE
*     IPHASE = 0
*     CALL NBLIST(I,RS2,NNB)
*     IPHASE = IP
*     DO 19 L = 2,NNB+1
*     JJ = ILIST(L)
*     LIST(L,2*IPAIR-1) = JJ
*  19 CONTINUE
*     LIST(1,2*IPAIR-1) = NNB
*
*       Save the perturbers for correction procedure.
      NNB = LIST(1,2*IPAIR-1)
      DO 20 L = 1,NNB
          J = LIST(L+1,2*IPAIR-1)
          JPERT(L) = J
   20 CONTINUE
*
*       Retain basic KS variables for explicit restart at merge termination.
      HM(IMERGE) = H(IPAIR)
      DO 25 K = 1,4
          UM(K,IMERGE) = U(K,IPAIR)
          UMDOT(K,IMERGE) = UDOT(K,IPAIR)
   25 CONTINUE
*
*       Save stellar type and binary flag.
      KSTARM(IMERGE) = KSTAR(I)
      IFLAGM(IMERGE) = LIST(2,2*IPAIR)
*
*       Set temporary KS components for hierarchy binary.
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
*
*       Obtain potential energy with respect to inner components.
      JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT1)
*
*       Save relative configuration.
      DO 30 K = 1,3
          XREL(K,IMERGE) = X(K,ICOMP) - X(K,JCOMP)
          VREL(K,IMERGE) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
*       Initialize primary velocity of JCOMP1 (needed in KSIN2).
          X0DOT(K,JCOMP1) = XDOT(K,JCOMP1)
   30 CONTINUE
*
*       Include interaction of inner c.m. & neighbours before making new KS.
      ICOMP = ICOMP1
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT2)
*
*       Copy mass of intruder to second KS component.
      BODY(JCOMP) = BODY(JCOMP1)
*
*       Replace name of #JCOMP with c.m. name (temporary for diagnostics).
      NAME2 = NAME(JCOMP)
      NAME(JCOMP) = NAME(JCOMP1)
*
*       Form new c.m. and initialize KS variables (JPERT safe first call!). 
      JCOMP = JCOMP1
      CALL KSIN2(1)
*
*       See whether modifications due to second binary are needed.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Initialize unperturbed ghost binary of outer component.
      T0(2*JPAIR-1) = 1.0D+06
      LIST(1,2*JPAIR-1) = 0
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
      CALL NBPOT(2,NNB,POT3)
      JLIST(1) = JCOMP1
      CALL NBPOT(1,NNB,POT4)
*
*       Update the merger energy to maintain conservation.
      EB1 = BODY(2*JPAIR-1)*BODY(2*JPAIR)*H(JPAIR)/BODY(JCOMP1)
      EMERGE = EMERGE + EB1
*
*       Save component masses and initialize ghost components.
      CM(3,IMERGE) = BODY(2*JPAIR-1)
      CM(4,IMERGE) = BODY(2*JPAIR)
      BODY(2*JPAIR-1) = 0.0D0
      BODY(2*JPAIR) = 0.0D0
*
*       Remove ghost from all perturber lists.
   50 JLIST(1) = JCOMP1
      ICM = N + KSPAIR
      CALL NBREM(ICM,1,NPAIRS)
*
*       Specify JCOMP1 as ghost of zero mass.
      BODY(JCOMP1) = 0.0D0
*
*       Initialize integration variables to prevent spurious predictions.
      DO 60 K = 1,3
          X0DOT(K,JCOMP1) = 0.0D0
          XDOT(K,JCOMP1) = 0.0D0
          F(K,JCOMP1) = 0.0D0
          FDOT(K,JCOMP1) = 0.0D0
          D2(K,JCOMP1) = 0.0D0
          D3(K,JCOMP1) = 0.0D0
   60 CONTINUE
*
*       Set large value of T0 which avoids integration of ghost particle.
      T0(JCOMP1) = 1.0D+06
*       Set large X0 & X to avoid perturber selection (no escape removal).
      X0(1,JCOMP1) = 1.0D+06
      X(1,JCOMP1) = 1.0D+06
*
*       Initialize c.m. & KS polynomials
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
      CALL KSIN2(2)
*
*       Define large negative c.m. name for identification & termination.
      NAME(ICM) = NAME(ICM) - 3*NZERO
*
*       Set c.m. & ghost names for merger identification (include escape).
      NAMEM(IMERGE) = NAME(ICM)
      NAMEG(IMERGE) = NAME(JCOMP1)
      NAME(JCOMP) = NAME2
*
*       Copy stability limit for termination test A(1 - E) < R0 in KSINT.
      R0(IPAIR) = PCRIT 
*
*       Update merger energy to maintain conservation.
      DPHI = (POT2 - POT1) + (POT4 - POT3)
      EMERGE = EMERGE + ZMU*HM(IMERGE) + DPHI
*
*       Specify non-zero indicator for sending new data to GRAPE in INTGRT.
      ISEND = -1
*       Set IPHASE < 0 for new sorting.
      IPHASE = -1
*
      RETURN
*
      END
