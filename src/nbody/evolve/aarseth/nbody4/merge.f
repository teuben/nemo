      SUBROUTINE MERGE
*
*
*       Merging of hierarchical binary.
*       -------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      INTEGER  JSAVE(5*LMAX)
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
*       Check merger limit.
      IF (NMERGE.GE.MMAX) THEN
          WRITE (6,5)  NMERGE
    5     FORMAT (/,5X,'WARNING!    MERGER LIMIT    NMERGE =',I4)
          GO TO 100
      END IF
*
*       Use double hierarchy procedure if one binary is already a merger.
      IF (NAME(N+KSPAIR).LT.0.OR.NAME(JCOMP).LT.0) THEN
          CALL MERGE2
          GO TO 100
      END IF
*
*       Check switching of inner and outer binary for improved restart.
      IF (JCOMP.GT.N) THEN
          SEMI2 = -0.5*BODY(JCOMP)/H(JCOMP-N)
*       Skip switch if #JCOMP is another merger.
          IF (SEMI2.GT.R(KSPAIR).AND.NAME(JCOMP).GT.0) THEN
              K = KSPAIR
              KSPAIR = JCOMP - N
              JCOMP = N + K
          END IF
      END IF
*
*       Set pair index & c.m. of inner binary and save outer component.
      IPAIR = KSPAIR
      I = N + IPAIR
      JCOMP1 = JCOMP
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
*       Save component masses and evaluate energy of the inner binary.
      CM(1,IMERGE) = BODY(2*IPAIR-1)
      CM(2,IMERGE) = BODY(2*IPAIR)
      CM(3,IMERGE) = 0.0D0
      ZMU = BODY(2*IPAIR-1)*BODY(2*IPAIR)/BODY(N+IPAIR)
*     EB = ZMU*H(IPAIR)
*
*       Form two-body elements of inner binary.
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      PMIN = SEMI*(1.0 - SQRT(ECC2))
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
*       Reduce index of outer component if JPAIR moved up at termination.
          IF (IPAIR.LT.JPAIR) THEN
              JPAIR = JPAIR - 1
              JCOMP1 = JCOMP1 - 1
          END IF
      END IF
*
*       Save the perturbers for correction procedure & rename if moved up.
      NNB = LIST(1,2*IPAIR-1)
      DO 20 L = 1,NNB
          J = LIST(L+1,2*IPAIR-1)
          IF (J.GT.I) J = J - 1
          IF (J.LE.2*NPAIRS.AND.J.GT.2*IPAIR) J = J - 2
          JSAVE(L) = J
*       Reduce steps of perturbers (new force uses c.m. approximation).
*         STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*       Include body #J in the time-step list unless already present.
*         IF (T0(J) + STEP(J).LT.TLIST) THEN
*             CALL NLMOD(J,1)
*         END IF
   20 CONTINUE
*
*       Delay saving KS variables in block version until end-point in KSTERM.
      IF (TIME.GT.TBLOCK) THEN
*       Retain basic KS variables for explicit restart at merge termination.
          HM(IMERGE) = H(IPAIR)
          DO 25 K = 1,4
              UM(K,IMERGE) = U(K,IPAIR)
              UMDOT(K,IMERGE) = UDOT(K,IPAIR)
   25     CONTINUE
      END IF
*
*       Save stellar type and binary flag.
      KSTARM(IMERGE) = KSTAR(I)
      IFLAGM(IMERGE) = LIST(2,2*IPAIR)
*
*       Terminate inner pair in order to merge components after updating.
      CALL KSTERM
*
*       Copy perturber list (NB! JPERT used by KSTERM).
      DO 28 L = 1,NNB
          JPERT(L) = JSAVE(L)
   28 CONTINUE
*
*       Obtain potential energy with respect to inner components.
      JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT1)
*
*       Save relative configuration and define old binary as composite body.
      DO 30 K = 1,3
          XREL(K,IMERGE) = X(K,ICOMP) - X(K,JCOMP)
          VREL(K,IMERGE) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
          X(K,ICOMP) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))
     &                                      /(BODY(ICOMP) + BODY(JCOMP))
          X0DOT(K,ICOMP) = (BODY(ICOMP)*X0DOT(K,ICOMP) +
     &                      BODY(JCOMP)*X0DOT(K,JCOMP))/
     &                                       (BODY(ICOMP) + BODY(JCOMP))
*       Initialize primary velocity of JCOMP1 (needed in KSREG & KSINIT).
          X0DOT(K,JCOMP1) = XDOT(K,JCOMP1)
   30 CONTINUE
*
*       Form new c.m. body & associated ghost of zero mass.
      BODY(ICOMP) = BODY(ICOMP) + BODY(JCOMP)
      BODY(JCOMP) = 0.0D0
*
*       Initialize integration variables to prevent spurious predictions.
      DO 40 K = 1,3
          X0DOT(K,JCOMP) = 0.0D0
          XDOT(K,JCOMP) = 0.0D0
          F(K,JCOMP) = 0.0D0
          FDOT(K,JCOMP) = 0.0D0
          D2(K,JCOMP) = 0.0D0
          D3(K,JCOMP) = 0.0D0
   40 CONTINUE
*
*       Set large value of T0 which avoids integration of ghost particle.
      T0(JCOMP) = 1.0D+06
*       Set large X0 & X to avoid perturber selection (no escape removal).
      X0(1,JCOMP) = 1.0D+06
      X(1,JCOMP) = 1.0D+06
*
*       See whether modifications due to second binary are needed.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Initialize unperturbed ghost binary of outer component.
      T0(2*JPAIR-1) = 1.0E+06
      LIST(1,2*JPAIR-1) = 0
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
*       Note that inner binary correction is included with IPAIR.
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
*       Include interaction of inner c.m. & neighbours to give net effect.
   50 JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT2)
*
*       Form square of c.m. velocity correction due to tidal effects.
      VI2 = X0DOT(1,ICOMP)**2 + X0DOT(2,ICOMP)**2 + X0DOT(3,ICOMP)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
      CORR = 1.0 + 2.0*DPHI/(BODY(ICOMP)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
*     DO 60 K = 1,3
*         X0DOT(K,ICOMP) = SQRT(CORR)*X0DOT(K,ICOMP)
*  60 CONTINUE
*
*       Remove ghost from all perturber lists.
      JLIST(1) = JCOMP
      CALL NBREM(ICOMP,1,NPAIRS)
*
*       Perform KS regularization of hierarchical system (ICOMP & JCOMP1).
      JCOMP = JCOMP1
      CALL KSREG
*
*       Define negative c.m. name for identification & termination.
      NAME(NTOT) = NZERO - NAME(NTOT)
*
*       Set c.m. & ghost names for merger identification (include escape).
      NAMEM(IMERGE) = NAME(NTOT)
      NAMEG(IMERGE) = NAME(JCOMP1)
*
*       Copy stability limit for termination test A(1 - E) < R0 in KSINT.
      R0(NPAIRS) = PCRIT 
*       Ensure that flag denotes new rather than primordial binary.
      LIST(2,JCOMP) = 0
*
*       Define merger energy to maintain conservation.
*     EMERGE = EMERGE + ZMU*HM(IMERGE)
*       Alternative: include tidal correction term and skip X0DOT adjustment.
      EMERGE = EMERGE + ZMU*HM(IMERGE) + DPHI
*
*       Set negative indicators for new sorting and sending data to GRAPE.
      IPHASE = -1
      ISEND = -1
*
  100 RETURN
*
      END
