      SUBROUTINE KSTERM
*
*
*       Termination of KS regularization.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CZ(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8  SAVE(14)
      CHARACTER*8  WHICH1
*
*
*       Copy pair index from COMMON save and define KS components & c.m.
      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
      JMIN = 0
      RIJ = R(IPAIR)
*       Reset indicator for enforcing new KS scheduling in routine SUBINT.
      NBPREV = 0
*
*       Prepare termination at block time (KS, triple, quad, chain or merge).
      IF (TIME.LE.TBLOCK.AND.IPHASE.LT.9) THEN
          TIME0 = TBLOCK
*       Skip KS integration for unperturbed orbit or T0(I1) at end of block.
          IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME0) GO TO 3
*
*       See whether the time interval should be modified by KSLOW procedure.
          IF (KSLOW(IPAIR).GT.1) THEN
              IMOD = KSLOW(IPAIR)
              ZMOD = FLOAT(ISLOW(IMOD))
          ELSE
              ZMOD = 1.0
          END IF
*
    1     DT0 = TIME0 - T0(I1)
*       Integrate up to current block time in case interval is too large.
          IF (DT0.GT.STEP(I1)) THEN
              TIME = T0(I1) + STEP(I1)
              H0(IPAIR) = H(IPAIR)
              Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
              CALL STUMPF(IPAIR,Z)
              CALL KSINT(I1)
              DTU = DTAU(IPAIR)
              STEP(I1) = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU
     &                                                   + R(IPAIR))*DTU
              STEP(I1) = ZMOD*STEP(I1)
*       Restrict increase of R for superfast particles in one block-step.
              IF (H(IPAIR).GT.100.0.AND.R(IPAIR).GT.RMIN) GO TO 3
              GO TO 1
          END IF
*
*       Determine the last regularized step by Newton-Raphson iteration.
          DTU = DT0/(R(IPAIR)*ZMOD)
          DTU = MIN(DTU,DTAU(IPAIR))
*       Include rare case of zero interval due to subtracting large values.
          DTU = MAX(DTU,1.0D-10)
          ITER = 0
    2     Y0 = DT0 - ZMOD*((ONE6*TDOT3(IPAIR)*DTU +
     &                             0.5*TDOT2(IPAIR))*DTU + R(IPAIR))*DTU
          YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
          YPR = ZMOD*YPR
          DTU = DTU - Y0/YPR
          DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
          DT1 = ZMOD*DT1
          ITER = ITER + 1
          IF (ABS(DT0 - DT1).GT.1.0E-10*STEP(I1).AND.ITER.LT.10) GO TO 2
*
*       Advance the KS solution to next block time and terminate at TIME0.
          DTAU(IPAIR) = DTU
          STEP(I1) = DT1
          TIME = T0(I1) + DT1
          H0(IPAIR) = H(IPAIR)
          Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
          CALL STUMPF(IPAIR,Z)
          CALL KSINT(I1)
    3     TIME = TIME0
*
*       Predict X & XDOT for body #JCOMP (note TIME = TBLOCK if second call).
          IF (JCOMP.GE.IFIRST) THEN
              CALL XVPRED(JCOMP,-1)
              IF (GAMMA(IPAIR).GT.0.2.AND.JCOMP.LE.N) THEN
                  JMIN = JCOMP
*       Initialize T0, X0 & X0DOT for XVPRED & FPOLY on large perturbation.
                  T0(JCOMP) = TIME
                  CALL DTCHCK(TIME,STEP(JCOMP),DTK(40))
                  DO 4 K = 1,3
                      X0(K,JCOMP) = X(K,JCOMP)
                      X0DOT(K,JCOMP) = XDOT(K,JCOMP)
    4             CONTINUE
              END IF
          END IF
      END IF
*
*       Predict coordinates and evaluate potential energy w.r.t. perturbers.
      CALL KSRES(IPAIR,J1,J2,0.0D0)
*
      NP = LIST(1,I1)
      DO 6 L = 1,NP
          JPERT(L) = LIST(L+1,I1)
    6 CONTINUE
*
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NP,POT1)
*
*       Rectify the orbit to yield U & UDOT consistent with binding energy.
      CALL KSRECT(IPAIR)
*
*       Retain final KS variables for explicit restart at merge termination.
      IF (TIME.LE.TBLOCK.AND.IPHASE.EQ.6) THEN
          HM(NMERGE) = H(IPAIR)
          RI = 0.0
          VI2 = 0.0
          DO 5 K = 1,4
              UM(K,NMERGE) = U(K,IPAIR)
              UMDOT(K,NMERGE) = UDOT(K,IPAIR)
              RI = RI + U(K,IPAIR)**2
              VI2 = VI2 + UDOT(K,IPAIR)**2
    5     CONTINUE
*
*       Include safety check in case of convergence problem at iteration.
          HI = (2.0*VI2 - (BODY(I1) + BODY(I2)))/RI
          IF (ABS(HI - H(IPAIR)).GT.0.01*ABS(HI)) THEN
              DM = BODY(I1) + BODY(I2) - BODY(ICM)
              WRITE (6,8)  IPAIR, NAME(I1), NAME(I2), R(IPAIR), RI,
     &                     H(IPAIR), HI, DM
    8         FORMAT (' DANGER!    KSTERM    KS NAM R U*U H0 HI DM ',
     &                                       I4,2I6,1P,5E10.2)
          END IF
      END IF
*
*       Check optional diagnostic output for disrupted new hard binary.
      IF (KZ(8).EQ.0) GO TO 10
      IF (LIST(2,I2).NE.0.OR.H(IPAIR).GT.0.0) GO TO 10
      IF (GAMMA(IPAIR).GT.0.5.AND.JCOMP.GT.0.OR.IPHASE.EQ.7) THEN
          IF (JCOMP.EQ.0.OR.IPHASE.EQ.7) JCOMP = I1
          K = 0
          IF (JCOMP.GT.N) THEN
              J2 = 2*(JCOMP - N)
              K = LIST(2,J2)
          END IF
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          EB = -0.5*BODY(I1)*BODY(I2)/SEMI
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          WRITE (8,9)  TTOT, NAME(I1), NAME(I2), K, NAME(JCOMP),
     &                 BODY(JCOMP), EB, SEMI, R(IPAIR), GAMMA(IPAIR), RI
    9     FORMAT (' END BINARY   T =',F6.1,'  NAME = ',2I5,I3,I6,
     &                      '  M(J) =',F8.4,'  EB =',F10.5,'  A =',F8.5,
     &                          '  R =',F8.5,'  G =',F5.2,'  RI =',F5.2)
      END IF
*
   10 IF (KZ(10).GT.1) THEN
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          WRITE (6,15)  TTOT, BODY(I1), BODY(I2), DTAU(IPAIR), R(IPAIR),
     &                  RI, H(IPAIR), IPAIR, GAMMA(IPAIR), STEP(I1),
     &                  LIST(1,I1)
   15     FORMAT (/,' END KSREG    TIME =',F8.2,2F8.4,F8.3,1PE10.1,
     &                                0PF7.2,F9.2,I5,F8.3,1PE10.1,0P,I5)
      END IF
*
*       Check optional output of degenerate binary (only new components).
      IF (KZ(8).GT.3.AND.(KSTAR(I1).GE.10.OR.KSTAR(I2).GE.10)) THEN
*       See whether either component has been regularized recently.
          NNB = LISTD(1) + 1
          K = 0
*       Check case of initial binary and loop over disrupted pairs.
          IF (IABS(NAME(I1) - NAME(I2)).EQ.1) K = -1
          DO 25 L = 2,NNB
              IF (NAME(I1).EQ.LISTD(L).OR.NAME(I2).EQ.LISTD(L)) K = -1
   25     CONTINUE
          IF (K.EQ.0.AND.H(IPAIR).LT.-0.1*ECLOSE) THEN
              CALL DEGEN(IPAIR,IPAIR,2)
          END IF
      END IF
*
*       Obtain global coordinates & velocities.
      CALL RESOLV(IPAIR,2)
*
*       Correct for differential potential energy due to rectification.
      CALL NBPOT(2,NP,POT2)
*       Add correction term with opposite sign for conservation.
*     ECOLL = ECOLL + (POT2 - POT1)
*     IF (ABS(POT1-POT2).GT.0.0001) WRITE (6,28)  POT1,BE(3),POT1-POT2
*  28 FORMAT (' CORRECT:    POT1 BE3 POT1-POT2  ',2F10.6,F10.6)
*
*       Reduce pair index, total number & single particle index.
      NPAIRS = NPAIRS - 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Save name of components & flag for modifying LISTD in UPDATE.
      JLIST(1) = NAME(I1)
      JLIST(2) = NAME(I2)
      JLIST(3) = LIST(2,I2)
*
*       Skip adjustment of tables if last or only pair being treated.
      IF (IPAIR.EQ.NPAIRS + 1) GO TO 60
*
*       Move the second component before the first.
      DO 50 KCOMP = 2,1,-1
          I = 2*IPAIR - 2 + KCOMP
*
          DO 30 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
   30     CONTINUE
*       Current velocity has been set in routine RESOLV.
          SAVE(7) = BODY(I)
          SAVE(8) = RADIUS(I)
          SAVE(9) = TEV(I)
          SAVE(10) = BODY0(I)
          SAVE(11) = EPOCH(I)
          SAVE(12) = TEV0(I)
          SAVE(13) = SPIN(I)
          SAVE(14) = ZLMSTY(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
          LAST = 2*NPAIRS - 1 + KCOMP
*
*       Move up global variables of other components.
          DO 40 J = I,LAST
              J1 = J + 1
              DO 35 K = 1,3
                  X(K,J) = X(K,J1)
   35         CONTINUE
              BODY(J) = BODY(J1)
              RADIUS(J) = RADIUS(J1)
              ZLMSTY(J) = ZLMSTY(J1)
              SPIN(J) = SPIN(J1)
              TEV(J) = TEV(J1)
              TEV0(J) = TEV0(J1)
              BODY0(J) = BODY0(J1)
              EPOCH(J) = EPOCH(J1)
              PHI(J) = PHI(J1)
              NAME(J) = NAME(J1)
              KSTAR(J) = KSTAR(J1)
              STEP(J) = STEP(J1)
              T0(J) = T0(J1)
   40     CONTINUE
*
*       Set new component index and copy basic variables.
          I = LAST + 1
          DO 45 K = 1,3
              X(K,I) = SAVE(K)
              X0DOT(K,I) = SAVE(K+3)
              XDOT(K,I) = SAVE(K+3)
   45     CONTINUE
          BODY(I) = SAVE(7)
          RADIUS(I) = SAVE(8)
          ZLMSTY(I) = SAVE(14)
          SPIN(I) = SAVE(13)
          TEV(I) = SAVE(9)
          TEV0(I) = SAVE(12)
          BODY0(I) = SAVE(10)
          EPOCH(I) = SAVE(11)
          NAME(I) = NAMEI
          KSTAR(I) = KSI
   50 CONTINUE
*
*       Update all regularized variables.
      CALL REMOVE(IPAIR,2)
*
*       Remove old c.m. from all COMMON tables (no F & FDOT correction).
      CALL REMOVE(ICM,3)
*
*       Set new global index of first & second component.
   60 ICOMP = 2*NPAIRS + 1
      JCOMP = ICOMP + 1
*
*       Modify all relevant COMMON list arrays.
      CALL UPDATE(IPAIR)
*
*       Form new force polynomial (skip triple, quad, merge & chain).
      IF (IPHASE.LT.4) THEN
          JFIRST = JCOMP + 1
*
*       Predict current coordinates & velocities for all other particles.
          CALL XVPRED(JFIRST,NTOT)
*
*       Obtain force polynomials (direct sum for critical cases).
          IF (RIJ.LT.0.01*RMIN) THEN
              CALL FPOLY1(ICOMP,JCOMP,2)
          ELSE
*       Obtain force polynomials on GRAPE (JCOMP is done if ICOMP = IFIRST).
              T0(JCOMP) = TIME
              CALL DTCHCK(TIME,STEP(JCOMP),DTK(40))
              CALL FPOLYI(ICOMP)
          END IF
*
*       Improve force polynomials of strong perturber after rectification.
          IF (JMIN.GE.IFIRST) THEN
              IPHASE = -1
              CALL FPOLYI(JMIN)
          END IF
*
*       See whether to include single components in NLIST.
          IF (T0(ICOMP) + STEP(ICOMP).LT.TLIST) THEN
              CALL NLMOD(ICOMP,1)
          END IF
          IF (T0(JCOMP) + STEP(JCOMP).LT.TLIST) THEN
              CALL NLMOD(JCOMP,1)
          END IF
      ELSE
*       Initialize T0 & X0 of old KS components to avoid prediction problem.
          T0(ICOMP) = TIME
          T0(JCOMP) = TIME
          CALL DTCHCK(TIME,STEP(ICOMP),DTK(40))
          CALL DTCHCK(TIME,STEP(JCOMP),DTK(40))
          DO 70 K = 1,3
              X0(K,ICOMP) = X(K,ICOMP)
              X0(K,JCOMP) = X(K,JCOMP)
   70     CONTINUE
      END IF
*
*       Check updating of neighbour list for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CMLIST
      END IF
*
      RETURN
*
      END
