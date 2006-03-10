
      SUBROUTINE MDOT
*
*
*       Mass loss from evolving stars.
*       ------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  M0,M1,LUM,MC,M10,MLWIND
      PARAMETER(MAXBLK=40)
      CHARACTER*8  WHICH1
      LOGICAL ICORR
      EXTERNAL MLWIND
*
*
*       Define current time (unit of million years).
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
      IPOLY = 0
      TIME0 = TIME
*
*       Find global index of next star to be checked (NB! reset KS & IKS).
    1 KS = 0
      IKS = 0
      IGHOST = 0
      IHDOT = 0
      IF (IPHASE.LT.0) IPOLY = -1
      IPHASE = 0
      TIME = TIME0
      DO 10 J = 1,NTOT
          IF (TEV(J).LE.TIME) THEN
              I = J
              GO TO 20
          END IF
   10 CONTINUE
*
*       Determine new evolution time (current TMDOT may be escaped star).
      KW = 0
      GO TO 70
*
*       Skip possible c.m. of compact subsystem or zero mass ghost.
   20 IF (NSUB.GT.0) THEN
          IF ((BODY(I).GT.BODY1.AND.NAME(I).GE.0).OR.NAME(I).EQ.0.OR.
     &        BODY(I).LE.0.0) THEN
              TEV(I) = TEV(I) + 0.01*TCR
              GO TO 1
          END IF
      END IF
*
*       Check for special treatment of possible hierarchical system (IHDOT).
      IF (I.LT.IFIRST) THEN
          KSPAIR = KVEC(I)
          ICM = N + KSPAIR
*       Select merger if #I is first or the second member is a c.m. body.
          IF (NAME(ICM).LT.0.AND.
     &       (I.EQ.2*KSPAIR - 1.OR.NAME(2*KSPAIR).GT.NZERO)) THEN
*       Activate indicator (distinguish inner and outer binary).
              IHDOT = 1
*       Exclude inner components and outer binary of double hierarchy.
              IF (NAME(ICM).LT.-2*NZERO) THEN
                  IF (I.LT.2*KSPAIR.OR.NAME(I).GT.NZERO) THEN
                      IPHASE = 7
                      CALL RESET
                      GO TO 1
                  END IF
              END IF
          END IF
      END IF
*
*       Determine relevant indices in case of merger.
      IF (IHDOT.GT.0) THEN
*
*       Identify ghost address and merger index from c.m.
          CALL FINDJ(I,IGHOST,IMERGE)
*
*       Determine index of mass-losing star for inner or outer binary.
          IF (IGHOST.LE.N) THEN
              IF (TEV(IGHOST).LT.TEV(I)) I = IGHOST
          ELSE
*       Save KS component for updating look-up time.
              I0 = I
*       Compare look-up times of #I and IGHOST <= N on first KS component.
              IF (I.LT.2*KSPAIR) THEN
                  IF (TEV(IGHOST).LT.TEV(I)) I = IGHOST
              ELSE
*       Form ghost pair index and select KS component by comparing TEV.
                  JPAIR = IGHOST - N
                  I = 2*JPAIR - 1
                  IF (TEV(2*JPAIR).LT.TEV(I)) I = 2*JPAIR
                  IHDOT = 2
                  WRITE (6,22)  IGHOST, I0, I, NAME(I), KSTAR(I), TEV(I)
   22             FORMAT (' OUTER GHOST:    IG I0 I NM K* TEV',5I5,F9.2)
                  IF (KSTAR(N+KSPAIR).GE.10) THEN
                      IPHASE = 7
                      CALL RESET
                      GO TO 1
                  END IF
              END IF
          END IF
      ELSE IF (BODY(I).LE.0.0D0.OR.NAME(I).EQ.0) THEN
*       Distinguish between merger ghost and member of CHAIN/TRIPLE/QUAD.
          IMERGE = 0
          ICM = 0
          J = I
*       Replace #I by corresponding ghost c.m. for zero-mass KS component.
          IF (I.LT.IFIRST) J = N + KVEC(I)
*
*       Identify merger index and determine corresponding c.m. body.
          DO 24 K = 1,NMERGE
              IF (NAMEG(K).EQ.NAME(J)) IMERGE = K
   24     CONTINUE
          IF (IMERGE.GT.0) THEN
              DO 25 K = N+1,NTOT
                  IF (NAMEM(IMERGE).EQ.NAME(K)) ICM = K
   25         CONTINUE
              KSPAIR = ICM - N
*       Terminate if ghost belongs to double hierarchy (bug fix 03/2001).
              IF (NAMEM(IMERGE).LT.-2*NZERO) THEN
                  IPHASE = 7
                  CALL RESET
                  GO TO 1
              END IF
          END IF
*
*       Define KS index and set indicator on successful identification.
          IF (ICM.GT.N) THEN
              IHDOT = 1
              IGHOST = J
              IF (IGHOST.GT.N) THEN
                  IHDOT = 2
                  JPAIR = IGHOST - N
                  I = 2*JPAIR - 1
                  IF (TEV(2*JPAIR).LT.TEV(I)) I = 2*JPAIR
                  I0 = IGHOST
              END IF
          ELSE
*       Extend TEV for other single ghosts or c.m. of compact subsystem.
              TEV(I) = TEV(I) + 0.01
              GO TO 1
          END IF
      END IF
*
*       Skip any c.m. particles after increasing look-up time.
      IF (I.GT.N) THEN
          IPAIR = I - N
          I1 = 2*IPAIR - 1
          I2 = I1 + 1
          TM = MIN(TEV(I1),TEV(I2))
          NWARN = NWARN + 1
          IF (NWARN.LT.100) THEN
              WRITE (6,30)  I, NAME(I), KSTAR(I1), KSTAR(I2), KSTAR(I),
     &                      TTOT, TM - TEV(I)
   30         FORMAT (' WARNING!    MDOT:    I NAM K* T D(TEV) ',
     &                              2I6,3I4,F9.2,1P,E10.2)
          END IF
          SEMI = -0.5*BODY(I)/H(IPAIR)
          IF (KZ(28).GT.0.AND.SEMI*SU.LT.10.0) THEN
              TEV(I) = TIME + 1.0/TSCALE
          ELSE
              TEV(I) = 1.0E+10
          END IF
          GO TO 1
      END IF
*
*       Set the initial mass and current type.
      NMDOT = NMDOT + 1
      M0 = BODY0(I)*ZMBAR
      KW = KSTAR(I)
*
*       Copy relevant mass (standard case or merger ghost member).
      IF (IGHOST.EQ.0.AND.IHDOT.EQ.0) THEN
          M1 = BODY(I)*ZMBAR
      ELSE
          IF (IHDOT.EQ.1) THEN
              K = 1
              IF (I.GE.IFIRST) K = 2
          ELSE
              K = 3
              IF (I.EQ.2*JPAIR) K = 4
          END IF
          M1 = CM(K,IMERGE)*ZMBAR
      END IF
      M10 = M1
*
*       Set interval for Reimers mass loss (small value on MS).
      DT = 1.0E+06*(TEV(I) - TEV0(I))*TSTAR
*
*       Obtain stellar parameters at previous epoch.
      AGE = TEV0(I)*TSTAR - EPOCH(I)
      MC = 0.0
      EPCH0 = EPOCH(I)
      CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC,RCC)
*
*       Ensure that type change occurs at time TEV.
      IF (KW.NE.KSTAR(I)) THEN
          KW = KSTAR(I)
          M0 = BODY0(I)*ZMBAR
          M1 = M10
          AGE = TEV0(I)*TSTAR - EPOCH(I)
          CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      END IF
*
*       Evaluate mass loss due to stellar wind during interval DT.
      IF (KW.LT.10) THEN
         DM = MLWIND(KW,LUM,RM,M1,MC,ZMET)*DT
      ELSE
         DM = 0.d0
      END IF
      DMR = ABS(DM/(M1 + 1.0d-10))
*
*       Restrict accumulated mass loss to maximum of 2 %.
      IF (DMR.GT.0.02) THEN
          DT = DT*0.02/DMR
          DM = 0.02*M1
          DMR = 0.02
      END IF
*
*       Check that mass loss does not exceed the envelope mass.
         IF (KSTAR(I).LT.10) THEN
            DML = MAX(M1 - MC,1.0D-07)
            IF (DML.LT.DM) THEN
               DT = (DML/DM)*DT
               DM = DML
*       Ensure correction for the last small mass loss (avoids looping).
               DMR = 0.011
            END IF
         END IF
*
*       Set indicator for mass loss correction (DM/M > 0.01).
      IF (DMR.GT.0.01) THEN
          ICORR = .TRUE.
          M1 = M1 - DM
*       Check rejuvenation of MS or HE star.
          IF (KW.LE.1.OR.KW.EQ.7) THEN
             XX = TM
             M0 = M1
             CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
             EPCH0 = TEV0(I)*TSTAR - (AGE*TM)/XX
          ELSE
             CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
          END IF
      ELSE
          ICORR = .FALSE.
      END IF
*
*       Set current age to actual look-up value (allows looping).
      TEVK = TEV0(I) + DT/(1.0D+06*TSTAR)
      AGE = TEVK*TSTAR - EPCH0
      AGE0 = AGE
*
      IF (AGE.GT.TN+0.01D0) THEN
*         WRITE (6,35)  NAME(I), KW, DM, AGE, AGE - TN
*  35     FORMAT (' WARNING!    AGE > TN    NM KW DMS AGE AGE-TN ',
*    &                                      I6,I4,F7.3,1P,2E10.2)
          IF (KW.LE.6) THEN
              AGE = 0.9999*TSCLS(11)
          ELSE
              AGE = 0.9999*TSCLS(5)
          END IF
      END IF
*
*       Obtain stellar parameters at current epoch.
      CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC,RCC)
*
*       Check inclusion of small accumulated mass loss on type change.
      IF (KW.NE.KSTAR(I)) THEN
          IF (KW.GE.10.OR.KW.EQ.7) THEN
              DM = M10 - M1
              ICORR = .TRUE.
          ELSE IF (.NOT.ICORR) THEN
              M1 = M1 - DM
              ICORR = .TRUE.
          END IF
      ELSE IF (KW.EQ.13.AND.I.LT.IFIRST) THEN
*       Ensure kick in FCORR for KS component unless large c.m. velocity.
          ICM =  N + KVEC(I)
          VI2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
*       Skip case of existing high velocity (i.e. TZ just after NS).
          IF (VI2.LT.100.0/VSTAR**2) THEN
               ICORR = .TRUE.
           END IF
      END IF
*
*       Set mass loss and new radius in N-body units.
      DMSUN = DM
      DM = DM/ZMBAR
      RNEW = RM/SU
      KW0 = KSTAR(I)
*
*       Check mass-loss treatment for inner binary components of merger.
      IF (IHDOT.GT.0) THEN
          IF ((KW.EQ.KW0.AND.AGE.LT.0.99*TN).OR.KW.LT.13) THEN
              IF (IHDOT.EQ.1) THEN
                  CALL HMDOT(IGHOST,IMERGE,M1,KW,MC,DMSUN,RNEW,ITERM)
              ELSE
                  CALL HMDOT2(IGHOST,IMERGE,M1,KW,MC,DMSUN,RNEW,ITERM)
              END IF
          ELSE
              ITERM = 1
          END IF
*
*       Terminate on non-zero indicator.
          IF (ITERM.GT.0) THEN
              IPHASE = 7
              CALL RESET
              GO TO 1
          END IF
      END IF
*
*       Include special procedures for KS components.
      IF (I.LT.IFIRST.AND.ICORR.AND.IHDOT.EQ.0) THEN
          SEMI = -0.5*BODY(N+KSPAIR)/H(KSPAIR)
*       Distinguish between KS pair and merger configuration.
          IF (NAME(N+KSPAIR).GT.0) THEN
*       Set random phase for neutron star or BH formation (negative index).
              IF (KW.EQ.13.OR.KW.EQ.14) THEN
                  JPAIR = -KSPAIR
                  KSTAR(I) = -KSTAR(I)
                  CALL KSAPO(JPAIR)
              END IF
*       Terminatenate for large mass loss or soft binary and re-determine index.
              IF ((DMR.GT.0.2.AND.R(KSPAIR).GT.RMIN).OR.
     &            (DM.GT.0.0.AND.H(KSPAIR) + DM/SEMI.GT.-ECLOSE).OR.
     &            KW.EQ.13.OR.KW.EQ.14) THEN
                  I = I + 2*(NPAIRS - KSPAIR)
*       Predict current KS variables and save at end of routine RESOLV.
                  CALL RESOLV(KSPAIR,3)
                  IPHASE = 2
                  JCOMP = 0
                  CALL KSTERM
                  KS = 1
              ELSE IF (DM.GT.0.0) THEN
*       Implement mass loss and expand KS orbit at constant eccentricity.
                  CALL HCORR(I,DM,RNEW)
              END IF
          ELSE
*       Adopt KS treatment for single outer component or terminate merger.
              IF (I.EQ.2*KSPAIR.AND.NAME(I).LE.NZERO.AND.
     &            NAME(N+KSPAIR).GT.-2*NZERO.AND.
     &            H(KSPAIR) + DM/SEMI.LT.-ECLOSE.AND.KW.NE.13) THEN
                  CALL HCORR(I,DM,RNEW)
              ELSE IF (DM.GT.0.0D0) THEN
                  IPHASE = 7
                  CALL RESET
                  GO TO 1
              END IF
          END IF
      END IF
*
*       Base new time scale for changes in radius & mass on stellar type.
      IF (KW.LE.1) THEN
         DTM = 0.05*TM
         DTR = TM - AGE
      ELSE IF (KW.EQ.2) THEN
         DTM = 0.05*(TSCLS(1) - TM)
         DTR = TSCLS(1) - AGE
      ELSE IF (KW.EQ.3) THEN
         IF (AGE.LT.TSCLS(6)) THEN
            DTM = 0.05*(TSCLS(4) - AGE)
         ELSE
            DTM = 0.05*(TSCLS(5) - AGE)
         END IF
         DTR = MIN(TSCLS(2),TN) - AGE
      ELSE IF (KW.EQ.4) THEN
         DTM = 0.05*TSCLS(3)
         DTR = MIN(TN,TSCLS(2) + TSCLS(3)) - AGE
      ELSE IF (KW.EQ.5) THEN
         IF (AGE.LT.TSCLS(9)) THEN
            DTM = 0.02*(TSCLS(7) - AGE)
         ELSE
            DTM = 0.02*(TSCLS(8) - AGE)
         END IF
         DTR = MIN(TN,TSCLS(13)) - AGE
      ELSE IF (KW.EQ.6) THEN
         IF (AGE.LT.TSCLS(12)) THEN
            DTM = 0.02*(TSCLS(10) - AGE)
         ELSE
            DTM = 0.02*(TSCLS(11) - AGE)
         END IF
         DTM = MIN(DTM,0.005d0)
         DTR = TN - AGE
      ELSE IF (KW.EQ.7) THEN
         DTM = 0.05*TM
         DTR = TM - AGE
      ELSE IF (KW.EQ.8.OR.KW.EQ.9) THEN
         IF (AGE.LT.TSCLS(6)) THEN
            DTM = 0.02*(TSCLS(4) - AGE)
         ELSE
            DTM = 0.02*(TSCLS(5) - AGE)
         END IF
         DTR = TN - AGE
      ELSE
*        DTM = AGE + 1.0
         DTM = 1.0E+04
         DTR = DTM
      END IF
*
*       Perform neighbour force corrections if mass loss is significant.
      IF (ICORR) THEN
*
*       Update the masses (single body, standard KS or merger).
          IF (IGHOST.EQ.0) THEN
              BODY(I) = M1/ZMBAR
          ELSE
              J = 2*KSPAIR - 2 + IHDOT
              BODY(J) = MAX(BODY(J) - DM,0.0D0)
          END IF
          BODY0(I) = M0/ZMBAR
*
*      Set new mass loss reference time TEV0.
          TEV0(I) = TEV(I)
*
*       Accumulate the total mass loss (solar units) and reduce cluster mass.
          ZMDOT = ZMDOT + DMSUN
          ZMASS = ZMASS - DM
          IF (KZ(19).GT.3) THEN
              WRITE (6,38)  I, NAME(I), KW, KSTAR(I), BODY(I)*ZMBAR,
     &                      DMSUN, ZMDOT, TPHYS
   38         FORMAT (' MDOT:   I NM KW K* MS DMS ZMDOT T6 ',
     &                          4I5,F6.1,F7.2,F7.1,F8.1)
          END IF
*
*       Replace any KS components by corresponding c.m. for main procedures.
          IF (I.LT.IFIRST) THEN
              IKS = I
              I = N + KSPAIR
              I1 = 2*KSPAIR - 1
              IF (LIST(1,I1).EQ.0) THEN
                  CALL RESOLV(KSPAIR,1)
              END IF
          END IF
*
*       Switch to the associated c.m. during force correction (save #I).
          IF (IGHOST.GT.0) THEN
              II = I
              I = N + KSPAIR
          END IF
*
*       Form neighbour list (depending on mass loss; IPHASE may be < 0).
          FACM = MIN(1.d0 + DMSUN,3.0d0)
          RS2 = (FACM*RSCALE)**2/FLOAT(N)**0.66667
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2
     &                                 + (X(3,I) - RDENS(3))**2
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
*       Adopt neighbour distance from interparticle spacing or density fit.
          RS2 = MAX(RS2,H2)
*       Restrict the guess for distant particles to avoid large NNB.
          IF (RI2.GT.4.0*RSCALE**2) THEN
              RI = SQRT(RI2)
              RS2 = (RI - 1.9*RSCALE)**2
              RS2 = MIN(RS2,0.01D0*RSCALE**2)
          END IF
          IF(IPHASE.NE.0) CALL GPSEND
          CALL NBLIST(I,RS2,NNB)
*       Ensure at least one neighbour.
          IF (NNB.EQ.0) THEN
              ILIST(2) = N
              IF (I.EQ.N) ILIST(2) = NTOT
              NNB = 1
          END IF
*       Include body #I at the end (counting from location #2).
          NNB2 = NNB + 2
          ILIST(NNB2) = I
*
*       Check prediction of neighbours and body #I to current time.
          DO 40 L = 2,NNB2
             J = ILIST(L)
             IF (T0(J).LT.TIME) THEN
                CALL XVPRED(J,-2)
             END IF
   40     CONTINUE
*
*       Specify IPHASE = -3 for using current members of ILIST in FPOLYI.
          IPHASE = -3
          ISEND = 0
          ILIST(1) = NNB
*
*       Perform total force & energy corrections (delay if DMSUN > 0.1).
          CALL FCORR(I,DM,KW)
*       Set flag to ensure new sorting after CALL FPOLYI (need IPHASE < 0).
          IF (I.GT.N) IPOLY = -1
*
*       Initialize new polynomials of neighbours & #I for DMSUN > 0.1.
          IF (ABS(DMSUN).GT.0.1.OR.KW.GE.13) THEN
*
*       Obtain new F & FDOT and time-steps (note FPOLYI may over-write I).
             DO 50 L = 2,NNB2
                J = ILIST(L)
                IF (L.EQ.NNB2) J = I
                CALL DTCHCK(TIME,STEP(J),DTK(MAXBLK))
                DO 45 KK = 1,3
                   X0DOT(KK,J) = XDOT(KK,J)
   45           CONTINUE
                CALL FPOLYI(J)
                IPOLY = -1
   50        CONTINUE
*       Use current potential for energy correction delayed from FCORR.
             ECDOT = ECDOT + DM*PHI(I)
             IF (KZ(14).GE.3.AND.KZ(31).GT.0) THEN
                CALL XVPRED(IFIRST,NTOT)
                CALL CMCORR
             END IF
             TPREV = TIME - STEPX
*       Set indicator < 0 for new sorting.
             IF (IGHOST.GT.0) I = II
             IF(IPHASE.EQ.-3) IPHASE = 0
          END IF
      END IF
*
*       Restore index in case of KS component (used as c.m. above).
      IF (IKS.GT.0) THEN
          I = IKS
*       Re-initialize KS polynomials for perturbed motion.
          IF (LIST(1,I1).GT.0) THEN
              CALL RESOLV(KSPAIR,1)
              CALL KSPOLY(KSPAIR,1)
*             TBLIST = TIME
          END IF
      END IF
*
*       Choose minimum of 0.05 of evolution time and remaining interval.
      DTM = MIN(DTM,DTR)
*
*       Impose a lower limit and convert time interval to scaled units.
      DTM = MAX(DTM,1.0D-04)/TSTAR
*
*       Update event counters & mass loss for new types > 2.
      if (kw.ne.kw0) then
         if (kw.eq.3) then
            nrg = nrg + 1
            zmrg = zmrg + dmsun
         else if (kw.eq.4) then
            nhe = nhe + 1
            zmhe = zmhe + dmsun
         else if (kw.eq.5) then
            nrs = nrs + 1
            zmrs = zmrs + dmsun
         else if (kw.ge.7.and.kw.le.9.and.kw0.le.6) then
            nnh = nnh + 1
            zmnh = zmnh + dmsun
         else if (kw.ge.10.and.kw.le.12.and.kw0.le.9) then
            nwd = nwd + 1
            zmwd = zmwd + dmsun
         else if ((kw.eq.13.or.kw.eq.15).and.kw0.le.12) then
            nsn = nsn + 1
            zmsn = zmsn + dmsun
         else if (kw.eq.14.and.kw0.le.13) then
            nbh = nbh + 1
            zmbh = zmbh + dmsun
         end if
      end if
*
*       Set new time for checking R & M and update R & classification type.
      TEV(I) = TEV(I) + DTM
      RADIUS(I) = RNEW
      KSTAR(I) = KW
*       Update epoch on transition to new type.
*     IF (KW.NE.KW0) THEN
          EPOCH(I) = AGE0 + EPCH0 - AGE
*     END IF
      IF (IHDOT.EQ.2) THEN
          TEV(I0) = TIME + DTM
*       Distinguish between single particle ghost and binary c.m.
          IF (NAME(IGHOST).LT.0.OR.NAME(IGHOST).GT.NZERO) THEN
              TEV(IGHOST) = TIME + DTM
          END IF
      END IF
*
*       Check optional diagnostic output.
      IF (KZ(19).GT.3.AND.(KW0.NE.KW.OR.DMR.GT.0.01)) THEN
          IF (KW0.NE.KW) THEN
              WHICH1 = ' TYPE   '
          ELSE
              WHICH1 = ' MASS   '
          END IF
          WRITE (6,65)  WHICH1, TPHYS, I, NAME(I), DMR, KW0, KW,
     &                  M0, M1, RADIUS(I)*SU, EMDOT
   65     FORMAT (' NEW',A8,' TPHYS I NAM DM/M KW0 KW M0 M R EMD ',
     &                        F7.1,2I5,F6.2,2I3,2F6.1,1P,E9.1,0P,F10.5)
      END IF
*
*       See if former KS pair can be regularized again.
      IF (KS.GT.0) THEN
          ICOMP = IFIRST
          JCOMP = IFIRST + 1
          RIJ2 = (X(1,ICOMP) - X(1,JCOMP))**2 +
     &           (X(2,ICOMP) - X(2,JCOMP))**2 +
     &           (X(3,ICOMP) - X(3,JCOMP))**2
          IF (RIJ2.LT.RMIN22) THEN
              CALL KSREG
*       Restore current time to prevent subsequent small steps.
              TIME = TBLOCK
              IF ((KW.EQ.13.OR.KW.EQ.14).AND.H(NPAIRS).LT.0.0) THEN
                  J = NTOT
                  SEMI = -0.5*BODY(J)/H(NPAIRS)
                  RA = R(NPAIRS)/SEMI
                  ECC2 = (1.0 - RA)**2 + TDOT2(NPAIRS)**2/(BODY(J)*SEMI)
                  WRITE (6,66)  KW, SQRT(ECC2), RA, SEMI*SU, STEP(J),
     &                          BODY(J)*ZMBAR, (XDOT(K,J)*VSTAR,K=1,3)
   66             FORMAT (' NS BINARY    KW E R/A A DT M V ',
     &                                   I4,2F6.2,1P,2E10.2,0P,4F7.1)
              END IF
          ELSE
              IPHASE = -1
          END IF
*       Ensure IPHASE < 0 on RETURN for new sorting after KSTERM.
          IPOLY = -1
          KS = 0
      END IF
*
*       Determine the time for next stellar evolution check.
   70 TMDOT = 1.0E+10
      DO 80 J = 1,NTOT
          IF (TEV(J).LE.TMDOT) THEN
              TMDOT = TEV(J)
          END IF
   80 CONTINUE
*
*       See if any other stars should be considered.
      IF (TMDOT.LT.TIME) GO TO 1
*
*       Update the maximum single body mass but skip compact subsystems.
      IF (NSUB.EQ.0) THEN
          BODY1 = 0.0
          DO 90 J = 1,N
              BODY1 = MAX(BODY1,BODY(J))
   90     CONTINUE
      END IF
*
*       Ensure IPHASE < 0 at the end and also enforce new KBLIST in SUBINT.
      IF (IPOLY.LT.0) THEN
          IPHASE = -1
          NBPREV = 0
      END IF
*
      RETURN
*
      END
