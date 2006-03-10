      SUBROUTINE IMPACT(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      CHARACTER*8  WHICH1
      REAL*8  XX(3,3),VV(3,3)
      INTEGER LISTQ(100)
      SAVE NMARG,LISTQ,QCHECK
      DATA NMARG,IZARE,LISTQ(1),QCHECK /0,0,0,0.0D0/
*
*
*       Set index of KS pair & both components of c.m. body #I.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      NTTRY = NTTRY + 1
      PERT1 = 0.0
      PERT2 = 0.0
      JCOMP = IFIRST
      NP = 0
      KS2 = 0
      RMAX2 = 1.0
      TTOT = TIME + TOFF
      IPX = 0
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
*
*       Copy perturber list (> 3 members) or obtain neighbour list.
      NNB2 = LIST(1,I1) + 1
      IF (NNB2.GT.4) THEN
          DO 2 L = 2,NNB2
              J = LIST(L,I1)
              NP = NP + 1
              JLIST(NP) = J
    2     CONTINUE
      ELSE
*       Restrict the initial guess for distant particles.
          IF (RI2.GT.4.0*RSCALE**2) THEN
              RI = SQRT(RI2)
              RS2 = (RI - 1.9*RSCALE)**2
              RS2 = MIN(RS2,0.01D0*RSCALE**2)
              RS2 = MAX(RS2,CMSEP2*R(IPAIR)**2)
*       Choose neighbour radius according to the perturbation.
          ELSE IF (GAMMA(IPAIR).LT.1.0E-04.AND.
     &             GAMMA(IPAIR).GT.1.0E-10) THEN
              RS2 = CMSEP2*MIN(R(IPAIR)**2,RMIN2)
              R2 = (2.0*BODYM/(BODY(I)*GAMMA(IPAIR)))**0.67*R(IPAIR)**2
              RS2 = MIN(R2,RS2)
              IF (LIST(1,I1).EQ.0) RS2 = 4.0*RS2
          ELSE IF (GAMMA(IPAIR).GT.0.1) THEN
              RS2 = 4.0*RMIN2
          ELSE IF (GAMMA(IPAIR).LT.1.0E-10) THEN
              RS2 = RSCALE**2/FLOAT(N)**0.66667
              H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
*       Adopt neighbour distance from interparticle spacing or density fit.
              RS2 = MAX(RS2,H2)
*       Note possible large RSCALE during final stages.
              R2 = 0.25*(RI2 + RC**2)
              RS2 = MIN(RS2,R2)
          ELSE
              RS2 = CMSEP2*MIN(R(IPAIR)**2,RMIN2)
          END IF
          ITER = 0
          IT1 = 0
*
*       Obtain the neighbour list from GRAPE (maximum 5 iterations).
    3     CALL NBLIST(I,RS2,NP)
*
*       Note counter for first member and repeat search if NP <= 1.
          IF (IT1.EQ.0.AND.NP.GT.0) IT1 = ITER
          IF (NP.LE.1) THEN
*       Add arbitrary perturber for probable hierarchy in the halo.
              IF (NP.GT.0.AND.RI2.GT.9.0*RSCALE**2) THEN
                  NP = NP + 1
                  IF (ILIST(2).LT.I-1) THEN
                      ILIST(NP+1) = ILIST(2) + 1
                  ELSE
                      ILIST(NP+1) = ILIST(2) - 1
                      IF (ILIST(NP+1).EQ.I) THEN
                          ILIST(NP+1) = NTOT
                      END IF
                  END IF
                  ITER = 5
              ELSE IF (ITER.GE.7) THEN
*       Increase gently after using fractional core radius.
                  RS2 = 2.0*RS2
                  ITER = ITER + 1
                  GO TO 3
              ELSE IF (NP.GT.0.AND.ITER - IT1.GE.2) THEN
                  ITER = 5
              END IF
              RS2 = 4.0*RS2
              ITER = ITER + 1
              IF (ITER.LT.4) GO TO 3
*       Try using the mean particle distance after four iterations.
              RS2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
              IF (ITER.LT.5) GO TO 3
*       Include one more search inside low-density core.
              IF (ITER.LE.6.AND.RI2.LT.RC**2) THEN
                  ITER = ITER + 1
                  RS2 = RC**2*MIN(10.0D0/FLOAT(NC),1.0D0)
                  GO TO 3
              ELSE IF (ITER.LE.6.AND.RI2.LT.RSCALE**2) THEN
                  ITER = ITER + 1
                  GO TO 3
              END IF
          END IF
*
*       Skip search if no perturber after four attempts.
          IF (NP.EQ.0.AND.ITER.GE.4) THEN
              IF (NDIAG.LT.9000) THEN
                  WRITE (18,4)  NAME(I), KSTAR(I), LIST(1,I1),
     &                          SQRT(RS2), R(IPAIR), SQRT(RI2),
     &                          GAMMA(IPAIR), STEP(I)
    4             FORMAT (' WARNING!    IMPACT    NM K* NP RS R r G DT',
     &                                            I7,2I4,1P,5E10.2)
                  CALL FLUSH(18)
              END IF
              NDIAG = NDIAG + 1
              GO TO 100
          END IF
*
*       Add an arbitrary perturber if NP = 1 after last call.
          IF (NP.EQ.1.AND.ITER.GE.4) THEN
              NP = NP + 1
              IF (ILIST(2).LT.I-1) THEN
                  ILIST(NP+1) = ILIST(2) + 1
              ELSE
                  ILIST(NP+1) = ILIST(2) - 1
                  IF (ILIST(NP+1).EQ.I) THEN
                      ILIST(NP+1) = NTOT
                  END IF
              END IF
              IPX = 1
          END IF
*
*       Copy members to JLIST.
          DO 5 L = 1,NP
              JLIST(L) = ILIST(L+1)
    5     CONTINUE
      END IF
*
*       Find the dominant body (JCOMP) and nearest perturber (JMAX).
      DO 10 L = 1,NP
          J = JLIST(L)
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          PERT = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PERT.GT.PERT2) THEN 
              IF (PERT.GT.PERT1) THEN
                  RJMIN2 = RIJ2
                  JMAX = JCOMP
                  JCOMP = J
                  PERT2 = PERT1
                  PERT1 = PERT
              ELSE
                  JMAX = J
                  PERT2 = PERT
                  RMAX2 = RIJ2
              END IF
          END IF
   10 CONTINUE
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Specify larger perturbation for optional chain regularization.
      IF (KZ(30).GT.0.AND.NCH.EQ.0) THEN
          GSTAR = 100.0*GMIN
          KCHAIN = 1
      ELSE
          GSTAR = GMIN
          KCHAIN = 0
      END IF
*
*       Only accept inward motion or small secondary perturbation.
      PERT3 = 2.0*R(IPAIR)**3*PERT2/BODY(I)
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 + 
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJ = SQRT(RJMIN2)
      A3 = 2.0/RIJ - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Set semi-major axis, eccentricity & apocentre of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      A0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      APO = ABS(SEMI)*(1.0 + ECC)
*
*       Quit on hyperbolic orbit with large impact parameter.
      IF (ECC1.GT.1.0.AND.PMIN.GT.50.0*SEMI) GO TO 100
*
*       Form binding energy of inner & outer binary.
      EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
      IF(ABS(EB).LT.1.0D-10) EB = -1.0D-10
      EB1 = -0.5*BODY(JCOMP)*BODY(I)/SEMI1
*
*       Obtain the total perturbing force acting on body #I & JCOMP.
      CALL FPERT(I,JCOMP,NP,PERT)
*
*       Choose maximum of dominant scalar & total vectorial perturbation.
      PERT = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      PERT4 = 2.0*RJMIN2*RIJ*PERT2/(BODY(I) + BODY(JCOMP))
      PERTM = MAX(PERT4,PERT)
*
*       Use combined semi-major axis for binary-binary collision.
      IF (JCOMP.GT.N) THEN
          JPAIR = JCOMP - N
          SEMI2 = -0.5D0*BODY(JCOMP)/H(JPAIR)
          J1 = 2*JPAIR - 1
          EB2 = -0.5*BODY(J1)*BODY(J1+1)/SEMI2
*       Define SEMI0 as smallest binary in case IPAIR denotes widest pair.
          SEMI0 = MIN(ABS(SEMI),ABS(SEMI2))
          SEMIX = MAX(SEMI,SEMI2)
          APO = APO + MAX(ABS(SEMI2),R(JPAIR))
          SEMI = SEMI + SEMI2
*       Do not allow negative or soft cross section.
          IF (1.0/SEMI.LT.0.5/RMIN) GO TO 100
*       Consider merger for PMIN > SEMI and large semi-major axis ratio.
          IF (PMIN.GT.SEMI.AND.SEMI2.GT.20.0*SEMI0) GO TO 30
      END IF
*
*       Check separation in case of chain regularization.
      IF (KCHAIN.GT.0) THEN
*       Form effective gravitational radius (combine triple & quad).
          EBT = EB + EB1
          ZMM = BODY(I1)*BODY(I2) + BODY(I)*BODY(JCOMP)
*       Set length of chain for decision-making (also used at termination).
          RSUM = R(IPAIR) + RIJ
          RI = R(IPAIR)
          IF (JCOMP.GT.N) THEN
              EBT = EBT + EB2
              ZMM = ZMM + BODY(J1)*BODY(J1+1)
              RSUM = RSUM + R(JPAIR)
              RI = MAX(R(JPAIR),RI)
          END IF
          RGRAV = ZMM/ABS(EBT)
*       Employ initial size as upper limit in case of weakly bound system.
          RGRAV = MIN(RGRAV,RMIN)
*       Save initial energy in binaries for routine SETSYS.
          EBCH0 = EBT - EB1
*       Use RIJ instead of RSUM in 3*RGRAV test (increases initial RIJ).
          IF (RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RSUM.GT.2.0*RMIN) GO TO 30
          GI = 2.0*BODY(JCOMP)*(RI/RIJ)**3/BODY(I)
*       Enforce KS orbit using MERGE for high eccentricity if PMIN > 10*RI.
          IF (ECC1.GT.0.99.AND.PMIN.GT.10.0*RI.AND.
     &        PERTM.LT.GMAX) GO TO 40
          IF (IPX.EQ.0.AND.GI.LT.0.02) GO TO 30
          IF (KZ(27).GT.0.AND.JCOMP.GT.N) THEN
              IF (SEMI0.LT.SEMI2) J1 = I1
              RT = 4.0*MAX(RADIUS(J1),RADIUS(J1+1))
*       Do not allow large distance ratio for nearly synchronous binary.
              IF (SEMI0.GT.RT.AND.RI.GT.25.0*SEMI0) GO TO 30
              IF (MIN(SEMI0,SEMI2).LT.0.05*RIJ) THEN
              IF (MAX(SEMI0,SEMI2).LT.0.1*RIJ) GO TO 30
              END IF
          END IF
      END IF
*
*       Include special case of strong interraction and large ECC1.
      IF (ECC1.GT.0.9.AND.GAMMA(IPAIR).GT.0.01) THEN
          IF (APO.LT.0.01*RMIN.AND.PMIN.LT.2.5*APO) GO TO 16
      END IF
*
*       Adopt triple, quad or chain regularization for strong interactions.
      IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.1.5*APO) GO TO 30
      IF ((RIJ.GT.RMIN.AND.SEMI1.GT.0.0).OR.RIJ.GT.2.0*RMIN) GO TO 100
      IF (PERTM.GT.100.0*GSTAR) GO TO 30
   16 IF (JCOMP.GT.N.AND.PMIN.GT.0.1*RMIN) THEN
          IF (PMIN.GT.A0 + SEMI2) GO TO 30
      END IF
      IF (JCOMP.GT.N.AND.PMIN.GT.4.0*SEMIX.AND.
     &   (ECC1.GT.0.9.AND.ECC1.LT.1.0)) GO TO 30
*
*       Check almost stable triples (factor 1.2 is experimental).
      IF (JCOMP.LE.N.AND.PMIN.GT.2.5*SEMI) THEN
          CALL HISTAB(IPAIR,JCOMP,PMIN,RSTAB)
          RA = SEMI1*(1.0 + ECC1)
          IF (SEMI1.LT.0.0) RA = RIJ
          GI = PERT*(RA/RIJ)**3
*       Use estimated apocentre perturbation for decision-making.
          IF (PMIN.GT.1.2*RSTAB) THEN
              IF (GI.LT.0.05) GO TO 30
*       Choose chain for critical case of highly eccentric outer orbit.
              IF (ECC1.LT.0.95) GO TO 100
          ELSE IF (PMIN.GT.0.9*RSTAB) THEN
*       Treat marginally stable triple according to external perturbation.
              IF (GI.LT.0.05) GO TO 30
              IF (GI.LT.1.0.OR.ECC1.LT.0.9) GO TO 100
          END IF
          IF (PMIN.GT.0.6*RSTAB.AND.PMIN.LT.0.9*RSTAB) GO TO 100
*       Delay for large distance ratio outside 0.5*RMIN.
          IF (RIJ.GT.MAX(10.0*APO,0.5*RMIN)) GO TO 100
          IF (RIJ.GT.10.0*APO) GO TO 100
          IF (PMIN.GT.3.0*APO) GO TO 100
      END IF
*
*       Skip chain if merged binary or chain c.m. (denoted by NAME <= 0).
      IF (NAME(I).LE.0.OR.NAME(JCOMP).LE.0) GO TO 100
*
*       Compare with existing subsystem of same type (if any).
      IF (NSUB.GT.0.AND.KCHAIN.EQ.0) THEN
          PERIM = R(IPAIR) + RIJ
          IF (JCOMP.GT.N) PERIM = PERIM + R(JPAIR)
          IGO = 0
          CALL PERMIT(PERIM,IGO)
          IF (IGO.GT.0) GO TO 100
      END IF
*
      WHICH1 = ' TRIPLE '
      IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
      IF (KCHAIN.GT.0) WHICH1 = ' CHAIN  '
*
      IF (H(IPAIR).GT.0.0) THEN
          WRITE (6,18)  I, JCOMP, ECC, ECC1, SEMI1, RIJ, GAMMA(IPAIR)
   18     FORMAT (' HYP CHAIN    I J E E1 A1 RIJ G  ',
     &                           2I6,2F7.3,1P,3E9.1)
      END IF
*
      IF (KZ(15).GT.1.OR.KZ(30).GT.1) THEN
          WRITE (6,20)  WHICH1, IPAIR, TTOT, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, LIST(1,I1)
   20     FORMAT (/,' NEW',A8,I4,'  T =',F8.2,'  H =',F6.0,
     &              '  R =',1P,E8.1,'  M =',0P,2F7.4,'  G4 =',1P,E8.1,
     &              '  R1 =',E8.1,'  P =',E8.1,'  E1 =',0P,F6.3,
     &              '  NP =',I2)
          CALL FLUSH(6)
      END IF
*
*       Include any close single or c.m. perturber (cf. routine SETSYS).
      IF (JMAX.NE.JCOMP.AND.SQRT(RMAX2).LT.MIN(2.0D0*RSUM,RMIN).AND.
     &    NAME(JMAX).GT.0) THEN
          IF (JCOMP.GT.N.AND.JMAX.GT.N) THEN
              JCMAX = 0
          ELSE
              WRITE (6,21)  NAME(JCOMP), NAME(JMAX), RSUM, SQRT(RMAX2)
   21         FORMAT (' B+2 CHAIN    NAM RSUM RMX ',2I7,1P,2E10.2)
              CALL XVPRED(JMAX,-1)
              JCMAX = JMAX
          END IF
      ELSE
          JCMAX = 0
      END IF
*
*       Save global index of intruder for TRIPLE or CHAIN.
      JCLOSE = JCOMP
*
*       Check B-B interaction for switch of IPAIR & JPAIR or inert binary.
      IF (KCHAIN.GT.0.AND.JCOMP.GT.N) THEN
          K1 = 2*JPAIR - 1
          WRITE (6,22)  NAME(I1), NAME(I2), NAME(K1), NAME(K1+1),
     &                  KSTAR(I), KSTAR(JCOMP), ECC, ECC1, A0, SEMI2,
     &                  RIJ, SEMI1, PMIN
   22     FORMAT (' CHAIN B-B    NAME K* E0 E1 A0 A2 RIJ A1 PM ',
     &                           4I6,2I4,2F7.3,1P,5E10.2)
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I2))
          IF (SEMI0.LT.4.0*RT.AND.LIST(1,J1).EQ.0.OR.
     &        MIN(SEMI0,SEMI2).LT.0.01*RIJ) THEN
*       Ensure that widest binary comes first (more similar to triple).
              IF (SEMI0.LT.SEMI2) THEN
                  KPAIR = JPAIR
                  JPAIR = IPAIR
                  IPAIR = KPAIR
                  JCLOSE = N + JPAIR
              END IF
*       Check reduction of c.m. index (JPAIR becomes JPAIR - 1 if > IPAIR).
              IF (JPAIR.GT.IPAIR) JCLOSE = JCLOSE - 1
              IF (KZ(26).LT.2) THEN
*       Replace unperturbed near-synchronous binary by inert body in CHAIN.
                  JCOMP = 0
                  WRITE (6,25)  SEMI0, RIJ, R(JPAIR), GAMMA(JPAIR)
   25             FORMAT (' INERT BINARY    A RIJ R G ',1P,4E10.2)
              END IF
          ELSE
              JCLOSE = 0
          END IF
      END IF
*
*       Set phase indicator for calling TRIPLE or QUAD from MAIN.
      IPHASE = 4
      KSPAIR = IPAIR
*
*       Include the case of two interacting KS pairs.
      IF (JCOMP.GT.N) THEN
          IPHASE = 5
*       Switch pair indices and rename JCOMP if JPAIR has smaller step.
          IF (STEP(J1).LT.STEP(I1).AND.LIST(1,I1).GT.0) THEN
              KSPAIR = JPAIR
              JCOMP = I
              KS2 = IPAIR
          ELSE
              KS2 = JPAIR
          END IF
          IF (KZ(27).LE.0.AND.JPAIR.GT.IPAIR) THEN
              IF (JCLOSE.GT.0) JCLOSE = JCLOSE - 1
          END IF
*       Reduce second index for later termination if higher.
          IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
      END IF
*
*       See whether chain regularization indicator should be switched on.
      IF (KCHAIN.GT.0) THEN
          IPHASE = 8
      END IF
*
*       Save KS indices and delay initialization until end of block step.
      CALL DELAY(KCHAIN,KS2)
*
*       Terminate binary in triple or widest binary-binary collision pair.
*     CALL KSTERM
*
*       Prepare procedure for chain between hierarchy and single body (9/99).
      IF (NAME(I).LT.0.AND.NAME(I).GE.-NZERO.AND.JCOMP.LE.N) THEN
*       Indentify merged ghost particle JG.
          CALL FINDJ(I1,JG,IM)
          WRITE (6,28)  NAME(I), NAME(JCOMP), NAME(JG),ECC1, PMIN, RIJ
   28     FORMAT (' HI CHAIN    NAM E1 PM RIJ ',I7,2I6,F7.3,1P,2E10.2)
          JJ = JCOMP
*       Terminate the merger in the usual way.
          KSPAIR = IPAIR
          IPHASE = 7
          CALL RESET
          ZMU = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)/BODY(NTOT)
          EBCH0 = EBCH0 + ZMU*H(NPAIRS)
*       Specify chain indicator and define the two single particles.
          IPHASE = 8
          JCMAX = JG
          JCLOSE = JJ
          KSPAIR = NPAIRS
*       Set relevant variables in DELAY before terminating inner binary.
          CALL DELAY(KCHAIN,KS2)
          CALL DELAY(IPHASE,-1)
*       Initialize new chain of the 4 members JMAX, JCLOSE & KS components.
          ISUB = 0
          CALL CHAIN(ISUB)
*       Note that IPHASE = -1 now and INTGRT goes back to the beginning.
      ELSE IF (NAME(I).LT.-NZERO.OR.NAME(JCOMP).LT.0.OR.
     &        (NAME(I).LT.0.AND.JCOMP.GT.N)) THEN
*       Continue until KS termination on MERGE2 or merger with JCOMP > N.
          IPHASE = 0
      END IF
*
      GO TO 100
*
*       Begin check for merger of stable hierarchical configuration.
   30 NMTRY = NMTRY + 1
      RA = SEMI1*(1.0 + ECC1)
      IF (SEMI1.LT.0.0) RA = RIJ
*
*       Identify formation of wide quadruple before merger is accepted.
      IF (JCOMP.GT.N.AND.ECC1.LT.1.0.AND.SEMI1.LT.0.1*RSCALE) THEN
          NNB = LISTQ(1) - 1
          K = 0
*       See whether current list contains first inner/outer component.
          NAM1 = NAME(2*JPAIR-1)
          DO 32 L = 2,NNB+2
              IF (NAM1.EQ.LISTQ(L)) K = K + 1
   32     CONTINUE
*       Generate diagnostics of first five outer orbits every half period.
          IF (K.LE.5.AND.TIME.GT.QCHECK.AND.KZ(18).GT.0) THEN
              ZMB = BODY(I) + BODY(JCOMP)
              RI = SQRT(RI2)
              TK = SEMI1*SQRT(SEMI1/ZMB)
              QCHECK = TIME + MIN(0.5*TWOPI*TK,0.1*TCR)
              TK = DAYS*TK
*       Employ the new stability criterion (MA 1997).
              Q = BODY(JCOMP)/BODY(I)
              XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
*       Include correction at small eccentricity (f(E) = 1.0 for now).
              FE  = 1.0
              PCR = 2.8*FE*XFAC**0.4*SEMIX
              WRITE (89,33)  TTOT, NAME(2*IPAIR-1), NAM1, K, RI,
     &                       ECC1, EB, EB2, EB1, TK, PMIN, PCR
   33         FORMAT (' QUAD    T NAM LQ RI E1 EB EB2 EB1 P1 PM PC ',
     &                          F8.1,2I6,I4,F6.2,F8.4,1P,3E12.3,3E9.1)
              CALL FLUSH(89)
*       Remove two oldest members if list is too big.
              IF (NNB.GT.96) THEN
                  DO 34 K = 2,NNB
                      LISTQ(K) = LISTQ(K+2)
   34             CONTINUE
                  NNB = NNB - 2
              END IF
*       Add current names (inner & outer) at end and update membership.
              LISTQ(NNB+3) = NAME(2*IPAIR-1)
              LISTQ(NNB+4) = NAME(2*JPAIR-1)
              LISTQ(1) = NNB + 3
          END IF
      END IF
*
*       Do not allow merger in the inner region of perturbed eccentric orbit.
      IF (RIJ.LT.SEMI1.AND.LIST(1,I1).GT.0) THEN
          IF (ECC1.GT.0.95.AND.RIJ.LT.2.0*PMIN) THEN
              GO TO 100
          END IF
      END IF
*
*       Allow temporary merger of inner part of extremely eccentric orbit.
      RFAC = 10.0*RMIN
      IF (ECC1.GT.0.99.AND.RA.GT.RFAC) THEN
          IF (RIJ.LT.0.1*SEMI1) RFAC = RA
      END IF
*
*       Increase apocentre tolerance to local scale factor for EB1 < EBS.
      EBS = 0.25*EBH/SQRT(1.0 + SQRT(RI2)/RSCALE)
      IF (EB1.LT.EBS) THEN
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
          RH = 6.0*SQRT(H2/CMSEP2)
          RFAC = MAX(RFAC,RH)
*       Extend maximum apocentre for massive systems (less perturbers).
          IF (BODY(I) + BODY(JCOMP).GT.10.0*BODYM) RFAC = 2.0*RFAC
      END IF
 
*       Skip merger for hyperbolic & soft binding energy or large apocentre.
      IF (EB.GT.EBH.OR.EB1.GT.EBS.OR.RA.GT.RFAC) THEN
          GO TO 100
      END IF
*
*       Estimate the relative apocentre perturbations on body #I & JCOMP.
      IF (ECC1.LT.0.95) THEN
          PERT = PERT*(RA/RIJ)**3
      ELSE
          PERT = PERT*(ABS(SEMI1)/RIJ)**3
      END IF
      PERTA = PERT4*(RA/RIJ)**3
*
*       Check tidal capture option (synchronous or evolving binary orbit).
      IF (KZ(27).GT.0) THEN
*       Skip merger if outer component would suffer tidal dissipation.
***       IF (SEMI1*(1.0 - ECC1).LT.4.0*RADIUS(JCOMP)) GO TO 100
*       Do not allow merger if Roche overflow or mass loss during next orbit.
          TK = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
          TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP),TEV(I))
*
*       Ensure SLEEP for circularizing binary with TCIRC > 1.0E+06.
          QPERI = A0*(1.0 - ECC)
          RM = MAX(RADIUS(I1),RADIUS(I2))
*
*       Delay merger for recently updated standard binary and short TCIRC.
          DT = MIN(TEV(I1),TEV(I2)) - TIME
          IF (KSTAR(I).EQ.0.AND.NAME(I).GT.0.AND.DT.LT.TK) THEN
              ICIRC = -1
              CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
*             WRITE (6,36) NAME(I1), ECC, TTOT, RADIUS(I1)*SU, QPERI, TC
*  36         FORMAT (' TCIRC    NAM E T R* QP TC ',
*    &                           I6,F7.3,F8.3,F7.1,1P,2E10.2)
*       Beware possible termination by routine HMDOT using QPERI < 3*RADIUS.
              IF (TC.LT.2000.0) GO TO 100
          END IF
          IF (KZ(19).GE.3) THEN
              IF (MIN(TEV(I1),TEV(I2)).LT.TIME + TK) GO TO 100
          END IF
*       Skip chaotic binary (KSTAR = -2 is treated below).
          IF (KSTAR(I).EQ.-1.OR.KSTAR(JCOMP).EQ.-1) GO TO 100
      END IF
*
*       Ensure consistency of estimated perturbations with termination.
      PERT = PERT + PERTA
      IF (NP.LE.3) THEN
          RS2 = CMSEP2*RMIN2
          CALL NBLIST(I,RS2,NP)
*       Copy neighbour list (routine FPERT skips body #JCOMP).
          DO 38 L = 1,NP
              JLIST(L) = ILIST(L+1)
   38     CONTINUE
      END IF
*
*       Evaluate the actual perturbation.
      CALL FPERT(I,JCOMP,NP,PERT2)
      PERT2 = PERT2*RJMIN2/(BODY(I) + BODY(JCOMP))
      GI = PERT2*(RA/RIJ)**3
*     IF (PERT4.GT.GMAX.OR.PERT.GT.0.1) GO TO 100
      IF (PERT4.GT.GMAX.OR.GI.GT.0.05) GO TO 100
 
*       Skip merger if an outer binary is fairly perturbed or not hard.
      IF (JCOMP.GT.N) THEN
          IF (GAMMA(JPAIR).GT.1.0E-03.OR.EB2.GT.EBH) GO TO 100
      END IF
*
*       Ensure the inner semi-major axis is used for subsequent tests.
   40 SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (RIJ.GT.5.0*RMIN) GO TO 100
*
*     -----------------------------------------------------------------------
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + ECC)*(1.0 + BODY(JCOMP)/BODY(I))**0.3333
*     FM = (2.0*BODY(JCOMP) - BODY(I))/(3.0*BODY(I))
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 82, 753) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*     -----------------------------------------------------------------------
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
*     QL = BODY(I)/BODY(JCOMP)
*     Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
*     Q3 = QL**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     EK = AR*SEMI*(1.0D0 + ECC)
*
*       Employ the new stability criterion (MA 1997).
      Q = BODY(JCOMP)/BODY(I)
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 1.0 + Q
      END IF
*       Include correction at small eccentricity (f(E) = 1.0 for now).
      FE  = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
*
*       Choose the most dominant triple in case of two binaries.
      YFAC = 1.0
      IF (JCOMP.GT.N) THEN
          SFAC = (1.0 + Q)**0.4*SEMI
          SFAC2 = (1.0 + BODY(I)/BODY(JCOMP))**0.4*SEMI2
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
          IF (SFAC2.GT.SFAC) THEN
              PCRIT = PCRIT*(SFAC2/SFAC)
          END IF
      END IF
*
*       Prepare inclination evaluation for triple or widest inner binary.
      IF (JCOMP.GT.N) THEN
*       Ensure widest inner binary (swap is OK for termination or ZARE).
          IF (SEMI.LT.SEMI2) THEN
              ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                             TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
              ECC = SQRT(ECC2)
              KPAIR = IPAIR
              IPAIR = JPAIR
              JPAIR = KPAIR
              I1 = 2*IPAIR - 1
              I2 = I1 + 1
              JJ = I
              I = JCOMP
              JCOMP = JJ
              SEMIZ = SEMI2
              SEMI2 = SEMI
              SEMI = SEMIZ
          END IF
      END IF
*
*       Resolve weakly perturbed binary (X(K,I1) = X(K,I2) for G < GMIN).
      IF (GAMMA(IPAIR).LT.GMIN) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Copy coordinates and velocities to local variables.
      DO 42 K = 1,3
          XX(K,1) = X(K,I1)
          XX(K,2) = X(K,I2)
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = XDOT(K,I1)
          VV(K,2) = XDOT(K,I2)
          VV(K,3) = XDOT(K,JCOMP)
  42  CONTINUE
*
*       Determine the inclination.
      CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
*
*       Adopt an emperical fudge factor for the inclination.
      IF (ECC1.LT.1.0) THEN
          YFAC = YFAC - 0.3*ANGLE/180.0
*       Employ additional gradual reduction above ECC1 = 0.96.
          IF (ECC1.GT.0.96) THEN
              YFAC = YFAC - 10.0*(ECC1 - 0.96)
*       Include enforcement for difficult case (experimental).
              IF (PMIN.LT.YFAC*PCRIT.AND.SEMI1.LT.0.2*RMIN) THEN
                  YFAC = 0.99*PMIN*(1.0 - PERT)/PCRIT
                  YFAC = MAX(YFAC,0.4D0)
              END IF
          END IF
      ELSE IF (RIJ.GT.RMIN) THEN
          GO TO 100
      END IF
*
*       Check whether the main perturber dominates the outer component.
      IF (JMAX.NE.JCOMP) THEN
          RIJ2 = (X(1,JMAX) - X(1,JCOMP))**2 +
     &           (X(2,JMAX) - X(2,JCOMP))**2 +
     &           (X(3,JMAX) - X(3,JCOMP))**2
          FMAX = (BODY(JMAX) + BODY(JCOMP))/RIJ2
          IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2) GO TO 100
      END IF
*
*       Include inclination angle procedure for marginal stability.
      IF (PMIN.LT.YFAC*PCRIT.AND.PMIN.GT.0.5*PCRIT.AND.
     &    LIST(1,I1).GT.0) THEN
          NMARG = NMARG + 1
*         YF = 1.0 - 0.3*ANGLE/180.0
          YF = 0.99*PMIN*(1.0 - PERT)/PCRIT
          TK = TWOPI*SEMI*SQRT(SEMI/(BODY(I) + BODY(JCOMP)))
*         WRITE (6,44)  NMARG, ANGLE, YF, PMIN, YF*PCRIT, TK
*  44     FORMAT (' MARGINAL    # ANGLE YF PM YF*PCR TK ',
*    &                          I7,2F7.2,1P,3E10.2)
*         CALL FLUSH(6)
          IF (PMIN*(1.0 - PERT).GT.YF*PCRIT.AND.
     &        (NMARG.GT.10000.OR.NMARG*TK.GT.0.5*DTADJ)) THEN
*    &        (NMARG.GT.100000.OR.NMARG*TK.GT.0.5*DTADJ)) THEN
              YFAC = 0.98*PMIN*(1.0 - PERT)/PCRIT
              WRITE (6,45)  NMARG, ANGLE, Q, YFAC, PMIN, PCRIT,
     &                      YFAC*PCRIT
   45         FORMAT (' NEW HIERARCHY    # IN Q YF PM PC1 PC2 ',
     &                                   I8,F7.1,2F6.2,1P,3E10.2)
              NMARG = 0
          END IF
      END IF
*
*       Determine time-scale for stability (absolute or approximate).
      PM1 = PMIN*(1.0 - 2.0*PERT)
      CALL TSTAB(I,ECC1,SEMI1,PM1,XFAC,YFAC,ITERM)
      IF (ITERM.GT.0) GO TO 100
*
*       Check perturbed stability condition (factor 1.01 avoids switching).
      IF (PMIN*(1.0 - PERT).LT.1.01*YFAC*PCRIT) GO TO 100
*
*       Check Zare exchange stability criterion and create diagnostics.
      IF (SEMI1.GT.0.0) THEN
          CALL ZARE(I1,I2,SP)
          IF (SP.LT.1.0.AND.ANGLE.LT.10.0) THEN
              IZARE = IZARE + 1
              IF (IZARE.LT.200) THEN
              WRITE (7,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
              WRITE (7,47) I,JCOMP,N,I1,I2,RIJ,SEMI1
   47         FORMAT (' I JCOMP N I1 I2 RIJ A1   ',5I6,1P,2E10.2)
              CALL FLUSH(7)
              WRITE (6,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
   48         FORMAT (' ZARE TEST    T Q E E1 A PM PCR YF SP ',
     &                               F8.2,F5.1,2F7.3,1P,3E9.1,0P,2F6.2)
              END IF
              GO TO 100
          END IF
          EK = 0.0
          WRITE (73,49)  TTOT, Q, ECC, ECC1, SEMI, PMIN, EK, PCRIT,
     &                   TG, SP, ANGLE, KSTAR(I)
   49     FORMAT (' STAB    T Q E E1 A PM EK PCR TG SP IN K* ',
     &                      F8.2,F5.1,2F7.3,1P,5E9.1,0P,F6.2,F7.1,I4)
          CALL FLUSH(73)
*         IF (KSTAR(I1).GE.10) TEV(I1) = 1.0E+10
*         IF (KSTAR(I2).GE.10) TEV(I2) = 1.0E+10
      END IF
*
*       Specify the final critical pericentre using the fudge factor.
      PCRIT = YFAC*PCRIT
*
      IF (NMERGE.EQ.MMAX) THEN
          IF (NWARN.LT.1000) THEN
              NWARN = NWARN + 1
              WRITE (6,50)  NMERGE
   50         FORMAT (5X,'WARNING!    MERGER LIMIT    NMERGE =',I4)
          END IF
*       Increase stability radii to enforce termination of weakest merger.
          DO 52 J = N+1,NTOT
              IF (NAME(J).LT.0) THEN
                  R0(J-N) = 1.1*R0(J-N)
              END IF
   52     CONTINUE
          GO TO 100
      END IF
*
*       Skip if #JCOMP is a chain c.m. but allow bound double hierarchy.
      IF (NAME(JCOMP).EQ.0) GO TO 100
      IF (ECC1.GT.1.0.AND.MIN(NAME(I),NAME(JCOMP)).LT.0) GO TO 100
*
*       Include diagnostics for double hierarchy or optional standard case.
      IF (NAME(I).LT.0.OR.NAME(JCOMP).LT.0) THEN
          WHICH1 = ' MERGE2 '
          WRITE (6,20)  WHICH1, IPAIR, TTOT, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, LIST(1,I1)
*       Note rare case of two hierarchies merging and identify ghost names.
           IF (NAME(I).LT.0.AND.NAME(JCOMP).LT.0) THEN
               CALL FINDJ(I1,JI,IM)
               CALL FINDJ(J1,JJ,JM)
               WRITE (6,60)  NAME(JI), NAME(JJ), ECC, ECC1, SEMI,
     &                       SEMI1, PMIN, PCRIT
   60          FORMAT (' HI MERGE    NAMG E E1 A A1 PM PC ',
     &                               2I6,2F7.3,1P,4E10.2)
           END IF
      ELSE IF (KZ(15).GT.1) THEN
          WHICH1 = ' MERGER '
          IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
          WRITE (6,20)  WHICH1, IPAIR, TTOT, H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, LIST(1,I1)
      END IF
*
*       Check for diagnostic output of quadruples.
      IF (SEMI1.GT.0.0.AND.JCOMP.GT.N.AND.KZ(18).GT.0) THEN
          ZMB = BODY(I) + BODY(JCOMP)
          TK = DAYS*SEMI1*SQRT(SEMI1/ZMB)
          WRITE (89,65)  TTOT, NAME(2*IPAIR-1), NAME(2*JPAIR-1),
     &                   LISTQ(1), SQRT(RI2), ECC1, EB, EB2, EB1,
     &                   TK, PMIN, PCRIT
   65     FORMAT (' QUAD#   T NAM LQ RI E1 EB EB2 EB1 P1 PM PC ',
     &                      F8.1,2I6,I4,F6.2,F8.4,1P,3E12.3,3E9.1)
      END IF
*
*       Generate a diagnostic file of stable hierarchies (suppressed).
      IF (ECC1.LT.-1.0) THEN
          RI = SQRT(RI2)/RC
          WRITE (80,70)  TPHYS, RI, NAME(JCOMP), QL, Q1, ECC, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, ANGLE, EMAX
   70     FORMAT (F8.1,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
          CALL FLUSH(80)
      END IF
*
*       Copy pair index and set indicator for calling MERGE from MAIN.
      KSPAIR = IPAIR
      IPHASE = 6
*
*       Save KS indices and delay merger until end of block step.
      CALL DELAY(KS2,KS2)
*
  100 IF (IPHASE.NE.8) JCMAX = 0
*
      RETURN
*
      END
