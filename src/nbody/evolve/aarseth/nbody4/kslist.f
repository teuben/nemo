      SUBROUTINE KSLIST(IPAIR)
*
*
*       KS perturber selection.
*       -----------------------
*
      INCLUDE 'common4.h'
*
*
*       Set component & c.m. index and form semi-major axis & eccentricity.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5d0*BODY(I)/H(IPAIR)
      EB = -0.5d0*BODY(I1)*BODY(I1+1)/SEMI
      ECC2 = (1.d0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2
     &                                   /(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Use semi-major axis and/or 4*RMIN for perturber selection.
      IF (EB.LT.EBH) THEN
          IF (IPHASE.EQ.6) THEN
              RAP = SEMI*(1.0 + ECC)
          ELSE
              RAP = SEMI*(1.d0 + ECC)
*       Allow increased search distance for energetic merged binaries.
              IF (NAME(I).LT.0) THEN
                  RAP = MAX(RAP,SEMI)
              END IF
          END IF
      ELSE
*       Include a tapered criterion depending on energy for soft binaries.
          IF (EB.LT.0.0) THEN
              ZFAC = 1.0 + ABS(EB - EBH)/ABS(EBH)
          ELSE
              ZFAC = 1.0
          END IF
*       Adopt actual apocentre for perturber selection if R > SEMI.
          IF (SEMI.GT.0.0.AND.R(IPAIR).GT.SEMI) THEN
              RAP = SEMI*(1.d0 + ECC)
          ELSE
              RAP = MAX(ZFAC*RMIN,R(IPAIR))
          END IF
*       Ensure extra perturbers at new regularizations (R may be small).
          IF (IPHASE.GT.0.AND.SEMI.GT.0.0D0) THEN
              RAP = MAX(RAP,SEMI)
              RAP = MIN(RAP,5.0D0*RMIN)
          END IF
      END IF
*
*       Assign search distance for neighbour list and tidal approximation.
      RAP2 = RAP**2
      RCRIT3 = 2.d0*RAP2*RAP/(BODY(I)*GMIN)
      FAC = MIN(1.d0 + BODY1/BODY(I),4.d0)
      RCRIT2 = FAC*CMSEP2*RAP2
*
*       Adopt cutoff from n(r) = N/2*(r/RSCALE)**2 to limit GRAPE list < LMAX.
      IF (IPHASE.NE.6.AND.N.GT.500) THEN
          RS2 = 2.d0*LMAX*RSCALE**2/FLOAT(N)
          RCRIT2 = MIN(RCRIT2,RS2)
      ELSE
          RCRIT2 = MIN(RCRIT2,RSCALE**2)
      END IF
*
*       Allow for wide binary perturber in differential correction.
      IF (RAP.LT.0.5*RMIN) THEN
          RCRIT2 = MAX(400.0*RMIN**2,RCRIT2)
      END IF
*
*       Impose minimum conditions for hyperbolic or soft binding energy.
      IF (EB.GT.0.0) THEN
          RCRIT2 = MIN(400.d0*RMIN2,0.09d0*RSCALE**2)
      ELSE IF (EB.GT.0.5*EBH) THEN
          RCRIT2 = MIN(RCRIT2,0.09*RSCALE**2)
      END IF
*
*       Determine neighbour list for c.m. body.
      CALL NBLIST(I,RCRIT2,NNB)
*
*       Select new perturbers (allow c.m. in KSPRED & 2 c.m. -> S in UPDATE).
      LM1 = LMAX - 3
    6 NNB1 = 1
      DO 10 L = 2,NNB+1
          J = ILIST(L)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
          RIJ3 = RSEP2*SQRT(RSEP2)
*       Estimate unperturbed distance from tidal limit approximation.
          IF (RIJ3.LT.BODY(J)*RCRIT3.AND.NNB1.LT.LM1) THEN
              NNB1 = NNB1 + 1
              LIST(NNB1,I1) = J
          ELSE IF (J.GT.N) THEN
*       Employ more generous criterion for possible wide binary.
              RJ = BODY(J)/H(J-N)
              IF (RSEP2.LT.CMSEP2*RJ**2.AND.NNB1.LT.LM1) THEN
                  NNB1 = NNB1 + 1
                  LIST(NNB1,I1) = J
              END IF
          END IF
   10 CONTINUE
*
*       Reduce perturber volume and repeat search if list is full.
      IF (NNB1.GE.LM1) THEN
          RCRIT3 = 0.9*RCRIT3*(FLOAT(LM1)/FLOAT(NNB))
          GO TO 6
      END IF
*
*       Ensure at least one perturber the first time.
      IF (NNB1.EQ.1.AND.IPHASE.GT.0.AND.NNB.GT.0) THEN
          RCRIT3 = 2.0*RCRIT3
          GO TO 6
      END IF
*
*       Check case of no perturbers (dual purpose).
      IF (NNB1.EQ.1) THEN
*       Add distant perturber for hyperbolic orbit (NB! check NNB = 0).
          IF (SEMI.LT.0.0) THEN
              NNB1 = 2
              LIST(2,I1) = ILIST(2)
              IF (NNB.EQ.0) THEN
                  LIST(2,I1) = IFIRST
              END IF
              GO TO 20
          END IF
*       Retain the previous perturber after partial reflection.
*         IF (TDOT2(IPAIR).GT.0.0D0.AND.KZ(25).GT.0) THEN
*             NNB1 = 2
*         ELSE IF (KZ(27).LE.0) THEN
          IF (KZ(27).LE.0) THEN
*       Specify one unperturbed period at apocentre (NB! check STEP(I)).
              STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
              STEP(I1) = MIN(STEP(I1),STEP(I))
          ELSE
*       Maintain perturbed motion during Chaos or Roche event.
              IF (KSTAR(I).EQ.-1) THEN
*             IF (KSTAR(I).EQ.-1.OR.(KSTAR(I).GT.10.AND.
*    &           MOD(KSTAR(I),2).EQ.1.AND.TEV(I).LT.TIME + STEP(I)))THEN
                  IF (LIST(1,I1).GT.0) THEN
                      NNB1 = 2
                      LIST(2,I1) = N
                  END IF
              ELSE
                  STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
                  STEP(I1) = MIN(STEP(I1),STEP(I))
              END IF
          END IF
*       Copy all neighbours (< LMAX-1) for soft binary in pericentre region.
          IF (R(IPAIR).LT.SEMI.AND.EB.GT.EBH) THEN
              DO 15 L = 2,NNB+1
                  IF (NNB1.LT.LM1) THEN
                      NNB1 = NNB1 + 1
                      LIST(NNB1,I1) = ILIST(L)
                  END IF
   15         CONTINUE
          END IF
      END IF
*
*       Save perturber membership.
   20 LIST(1,I1) = NNB1 - 1
*
*     IF (IPHASE.EQ.6)  WRITE (6,30) NAME(I1), NNB, NNB1-1, RAP
*  30 FORMAT (' KSLIST (MERGER)    NM NNB1 NP RAP  ',I6,2I4,1P,E10.2)
      RETURN
*
      END
