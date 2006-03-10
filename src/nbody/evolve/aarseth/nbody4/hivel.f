      SUBROUTINE HIVEL(IH)
*
*
*       High-velocity particle search.
*       ------------------------------
*
      INCLUDE 'common4.h'
      SAVE NOUT
      DATA NOUT /0/
*
*
*       Check for removal of distant high-velocity particle.
      IF (IH.LT.0) THEN
    1     LH = 0
          DO 2 L = 1,NHI
              I = LISTV(L)
              VI2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
              RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                       (X(3,I) - RDENS(3))**2
*       Save index of fast particle outside 3*<R>.
              IF (RI2.GT.9.0*RSCALE**2) THEN
                  LH = L
                  LI = I
                  RL2 = RI2
                  VL2 = VI2
              END IF
    2     CONTINUE
*       Reduce membership and remove any distant member from the list.
          IF (LH.GT.0) THEN
              WRITE (29,3)  LI, NAME(LI), SQRT(RL2), SQRT(VL2)
    3         FORMAT (' HIVEL REMOVE    I NAM R V ',2I6,1P,2E10.1)
              CALL FLUSH(29)
              NHI = NHI - 1
              DO 4 L = LH,NHI
                  LISTV(L) = LISTV(L+1)
    4         CONTINUE
              GO TO 1
          END IF
      END IF
*
*       Specify square velocity limit in terms of current state.
      VMAX2 = 16.0*ECLOSE
*
*       Set index to known neutron star or terminated KS/chain.
      IF (IH.GT.0) THEN
          I1 = IH
          I2 = IH
      ELSE
*       Include two first single particles or up to five chain members.
          I1 = IFIRST
          I2 = IFIRST + 1
          IF (IPHASE.EQ.8) THEN
              I2 = IFIRST + 4
*       Search all particles after escaper removal.
          ELSE IF (IH.EQ.0.AND.IPHASE.EQ.-2) THEN
              I2 = NTOT
              NHI = 0
          ELSE IF (IPHASE.EQ.1) THEN
*       See whether the first few locations may have been exchanged.
              DO 5 L = 1,NHI
                  IF (LISTV(L).LT.IFIRST + 3) THEN
                      I2 = NTOT
                  END IF
    5         CONTINUE
          END IF
      END IF
*
*       Add any new high-velocity particles (save F**2 > N & STEP < DTMIN).
      NHV = 0
      DO 10 I = I1,I2
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
          VI2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
*       Form a list of recently ejected candidates.
          IF (FI2.GT.FLOAT(N).OR.STEP(I).LT.DTMIN) THEN
              IF (VI2.GT.VMAX2) THEN
                  NHV = NHV + 1
                  JLIST(NHV) = I
              END IF
              GO TO 10
          END IF
          IF (VI2.GT.VMAX2.AND.RI2.LT.9.0*RSCALE**2) THEN
              DO 8 L = 1,NHI
                  IF (I.EQ.LISTV(L)) GO TO 10
    8         CONTINUE
*       Check maximum membership and possible ghost particle.
              IF (NHI.GE.20.OR.STEP(I).GT.1.0) GO TO 10
              NHI = NHI + 1
              LISTV(NHI) = I
              NOUT = NOUT + 1
              WRITE (29,9)  TTOT, NHI, I, NAME(I), KSTAR(I), SQRT(VI2),
     &                                                SQRT(RI2), STEP(I)
    9         FORMAT (' HIVEL ADD    T NHI I NM K* VI R DT ',
     &                               F10.4,I4,2I6,I4,2F6.2,1P,E10.2)
              CALL FLUSH(29)
          END IF
   10 CONTINUE
*
*       Consider single fast particle or hyperbolic two-body motion.
      IF (NHI.LT.19.AND.(NHV.EQ.1.OR.NHV.EQ.2)) THEN
*       Compare any candidates with existing members.
          DO 20 K = 1,NHV
              DO 15 L = 1,NHI
                  IF (JLIST(K).EQ.LISTV(L)) GO TO 30
   15         CONTINUE
   20     CONTINUE
          I1 = JLIST(1)
*       Include single particles without further tests.
          IF (NHV.EQ.1) THEN
              NHI = NHI + 1
              LISTV(NHI) = I1
              WRITE (29,22)  TTOT, NHI, NAME(I1), IPHASE, STEP(I1)
   22         FORMAT (' HIVEL ADD    T NHI NM IPH DT ',
     &                               F10.4,I4,I6,I4,1P,E10.2)
              GO TO 30
          END IF
*       Evaluate two-body energy.
          I2 = JLIST(2)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 25 K = 1,3
              RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
              VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
              RDOT = RDOT + (X(K,I1) - X(K,I2))*(XDOT(K,I1)-XDOT(K,I2))
   25     CONTINUE
          RIJ = SQRT(RIJ2)
          SEMI = 2.0/RIJ - VIJ2/(BODY(I1) + BODY(I2))
*       Accept outwards hyperbolic motion arising from recent interaction.
          IF (SEMI.LT.0.0.AND.RIJ.LT.10.0*RMIN.AND.RDOT.GT.0.0) THEN
              NHI = NHI + 1
              LISTV(NHI) = I1
              NHI = NHI + 1
              LISTV(NHI) = I2
              WRITE (29,28)  TTOT, NHI, NAME(I1), NAME(I2), IPHASE, RIJ
   28         FORMAT (' HIVEL ADD    T NHI NM IPH RIJ',
     &                               F10.4,I4,2I6,I4,1P,E10.2)
          END IF
      END IF
*
*       Increase counter if loop is over all particles.
   30 IF (I2.EQ.NTOT) THEN
          NHIVEL = NHIVEL + 1
      END IF
*
      RETURN
*
      END
