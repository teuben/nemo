      SUBROUTINE KSINIT
*
*
*       Initialization of KS regularization.
*       ------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  Q(3),RDOT(3),UI(4),VI(4),A1(3,4)
*
*
*       Set new global indices of the components and current pair index.
      ICOMP = 2*NPAIRS - 1
      JCOMP = ICOMP + 1
      IPAIR = NPAIRS
*
*       Specify mass & name for new c.m. and initialize radius & type.
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
      NAME(NTOT) = NZERO + NAME(ICOMP)
      RADIUS(NTOT) = 0.D0
      ZLMSTY(NTOT) = 0.D0
      SPIN(NTOT) = 0.D0
      TEV(NTOT) = 1.0E+10
      TEV0(NTOT) = TIME
      BODY0(NTOT) = BODYM
      EPOCH(NTOT) = TIME*TSTAR
      KSTAR(NTOT) = 0
      IMOD = 1
*
*       Ensure that both components are updated at the same time.
*      IF(NAME(ICOMP).LE.N.AND.NAME(JCOMP).LE.N)THEN
*         TEV(ICOMP) = MIN(TEV(ICOMP),TEV(JCOMP))
*         TEV(JCOMP) = TEV(ICOMP)
*      ENDIF
*
*       Define c.m. coordinates & velocities and set XDOT for components.
      DO 10 K = 1,3
          X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(NTOT)
          X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(NTOT)
          XDOT(K,NTOT) = X0DOT(K,NTOT)
          XDOT(K,ICOMP) = X0DOT(K,ICOMP)
          XDOT(K,JCOMP) = X0DOT(K,JCOMP)
   10 CONTINUE
*
*       Skip KS initialization at merger termination (H, U & UDOT in RESET).
      IF (IPHASE.EQ.7) THEN
          EB = 2.0*EBH
          GO TO 50
      END IF
*
*       Define relative coordinates and velocities in physical units.
      DO 20 K = 1,3
          Q(K) = X(K,ICOMP) - X(K,JCOMP)
          RDOT(K) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
   20 CONTINUE
*
*       Introduce regularized variables using definition of 1985 paper.
      R(IPAIR) = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2)
*
*       Initialize the regularized coordinates according to sign of Q(1).
      IF (Q(1).LE.0.0D0) THEN
          UI(3) = 0.0D0
          UI(2) = SQRT(0.5D0*(R(IPAIR) - Q(1)))
          UI(1) = 0.5D0*Q(2)/UI(2)
          UI(4) = 0.5D0*Q(3)/UI(2)
      ELSE
          UI(4) = 0.0D0
          UI(1) = SQRT(0.5D0*(R(IPAIR) + Q(1)))
          UI(2) = 0.5D0*Q(2)/UI(1)
          UI(3) = 0.5D0*Q(3)/UI(1)
      END IF
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Form regularized velocity and set initial KS coordinates & TDOT2.
      TDOT2(IPAIR) = 0.0D0
      DO 30 K = 1,4
          UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &                                                  A1(3,K)*RDOT(3))
*       Note that A1(J,K) is the transpose of A1(K,J).
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = U(K,IPAIR)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*UI(K)*UDOT(K,IPAIR)
   30 CONTINUE
*
*       Evaluate initial binding energy per unit mass and EB.
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                   UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                                              BODY(NTOT))/R(IPAIR)
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
*
*       Obtain force polynomial for c.m. with components ICOMP & JCOMP.
   50 IF (R(IPAIR).LT.0.01*RMIN) THEN
          CALL FPOLY1(ICOMP,JCOMP,1)
*       NB! Needs GPSEND at TIME = 0 (bug 14-6-02).
          CALL GPSEND
      ELSE
          CALL FPOLYI(NTOT)
      END IF
*
*       See whether to include new c.m. or single components in NLIST.
      IF (T0(NTOT) + STEP(NTOT).LT.TLIST) THEN
          CALL NLMOD(NTOT,1)
      END IF
*
*       Form perturber list.
      CALL KSLIST(IPAIR)
*
*       Transform any unperturbed hard binary to apocentre and set time-step.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      IF (LIST(1,ICOMP).EQ.0.AND.EB.LT.EBH) THEN
          TK = TWOPI*ABS(SEMI)*SQRT(ABS(SEMI)/BODY(NTOT))
*       Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
          IF (IPHASE.NE.7.AND.IPHASE.NE.8) THEN
              DO 55 K = 1,4
                  VI(K) = UDOT(K,IPAIR)
   55         CONTINUE
*       Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
              CALL TPERI(SEMI,UI,VI,BODY(NTOT),TP)
              STEP(ICOMP) = 0.5*MIN(TK,STEP(NTOT)) - TP
*       Transform KS variables to peri and by pi/2 to apocentre (skip apo).
              IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
                  TIME0 = TIME
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
                  TIME = TIME0
              ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
                  TDOT2(IPAIR) = -1.0E-20
              END IF
          END IF
      END IF
*
*       Estimate an appropriate KS slow-down index for G < GMIN.
      IF (LIST(1,ICOMP).EQ.0.AND.SEMI.GT.0) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(NTOT))
          IF (KZ(26).GT.0.AND.STEP(NTOT).GT.TK) THEN
              IMOD = 1 + LOG(STEP(NTOT)/TK)/0.69
              IMOD = MIN(IMOD,5)
          END IF
      END IF
*
*       Specify zero membership and large step for second component.
      LIST(1,JCOMP) = 0
      STEP(JCOMP) = 1.0E+06
*
*       Obtain polynomials for perturbed KS motion (standard case).
      CALL KSPOLY(IPAIR,IMOD)
*
*       Set maximum of RMIN & 2*SEMI as termination scale for hard binary.
      IF (EB.LT.EBH) THEN 
          R0(IPAIR) = MAX(RMIN,2.0D0*SEMI)
          IF (R(IPAIR).LT.0.25D0*SEMI) THEN
              R0(IPAIR) = MIN(R0(IPAIR),4.D0*RMIN)
          END IF
      ELSE
          R0(IPAIR) = R(IPAIR)
      END IF
*
*       Increase regularization counters (NKSHYP for hyperbolic orbits).
      NKSREG = NKSREG + 1
      IF (H(IPAIR).GT.0.0) NKSHYP = NKSHYP + 1
*
      IF (KZ(10).GT.0) THEN
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          WRITE (6,60)  TTOT, NAME(ICOMP), NAME(JCOMP), DTAU(IPAIR),
     &                  R(IPAIR), RI, H(IPAIR), IPAIR, GAMMA(IPAIR),
     &                  STEP(NTOT), LIST(1,ICOMP)
   60     FORMAT (/,' NEW KSREG    TIME =',F8.2,2I6,F12.3,1P,E10.1,0P,
     &                             F7.2,F9.2,I5,F8.3,1P,E10.1,0P,I5)
      END IF
*
*       Modify the termination criterion according to value of NPAIRS.
      IF (NPAIRS.GT.KMAX - 3) GMAX = 0.8*GMAX
      IF (NPAIRS.LT.KMAX - 5.AND.GMAX.LT.0.001) GMAX = 1.2*GMAX
      IF (NPAIRS.EQ.KMAX) WRITE (6,70)  NPAIRS, TTOT
   70 FORMAT (5X,'WARNING!    MAXIMUM KS    NPAIRS TIME',I5,F8.2)
*
*       See whether either component has been regularized recently.
      NNB = LISTD(1) + 1
      K = 0
*       Check case of initial binary and loop over disrupted pairs.
      IF (IABS(NAME(ICOMP) - NAME(JCOMP)).EQ.1) K = -1
      DO 80 L = 2,NNB
          IF (NAME(ICOMP).EQ.LISTD(L).OR.NAME(JCOMP).EQ.LISTD(L)) K = -1
   80 CONTINUE
*
*       Check optional output of degenerate binary (only new components).
      IF (KZ(8).GT.3.AND.MAX(KSTAR(ICOMP),KSTAR(JCOMP)).GE.10) THEN
            IF (K.EQ.0.AND.H(IPAIR).LT.-0.1*ECLOSE.AND.IPHASE.NE.6) THEN
                CALL DEGEN(IPAIR,IPAIR,1)
          END IF
      END IF
*
*       Ensure that mergers are treated as new binaries.
      IF (IPHASE.EQ.6) K = 0
*       Set flags to distinguish primordial binaries & standard KS motion.
      LIST(2,JCOMP) = K
      KSLOW(IPAIR) = 1
*
*       Check diagnostic output of new hard binary.
      IF (KZ(8).GT.0.AND.K.EQ.0) THEN
          IF (EB.GT.EBH) GO TO 90
          SEMI = -0.5*BODY(NTOT)/H(IPAIR)
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          IF (IPHASE.EQ.6) K = -1
          WRITE (8,85)  TTOT, NAME(ICOMP), NAME(JCOMP), K, BODY(ICOMP),
     &                  BODY(JCOMP), EB, SEMI, R(IPAIR), GAMMA(IPAIR),
     &                  RI
   85     FORMAT (' NEW BINARY   T =',F7.1,'  NAME = ',2I5,I3,
     &                        '  M =',2F8.4,'  EB =',F9.4,'  A =',F8.5,
     &                          '  R =',F8.5,'  G =',F6.3,'  RI =',F5.2)
      END IF
*
*       Include diagnostics for new hard binary and degenerate component(s).
   90 IF (MAX(KSTAR(ICOMP),KSTAR(JCOMP)).GE.10.AND.K.EQ.0.AND.
     &    H(IPAIR).LT.-ECLOSE.AND.NAME(JCOMP).LE.NZERO) THEN
          SEMI = -0.5*BODY(NTOT)/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                 TDOT2(IPAIR)**2/(BODY(NTOT)*SEMI)
          WRITE (6,95)  TTOT, IPAIR, NAME(ICOMP), NAME(JCOMP),
     &                  KSTAR(ICOMP), KSTAR(JCOMP), SQRT(ECC2), SEMI,
     &                  GAMMA(IPAIR), LIST(1,ICOMP)
   95     FORMAT (' NEW DEGEN    T KS NAM K* E A G NP ',
     &                           F9.2,I5,2I6,2I4,F7.3,1P,2E10.2,0P,I4)
      END IF
*
      RETURN
*
      END
