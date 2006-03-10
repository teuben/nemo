      SUBROUTINE KSTIDE(IPAIR,QPERI)
*
*
*       Tidal interaction of KS pair.
*       -----------------------------
*
      INCLUDE 'common4.h'
      REAL*8  DE(2)
      CHARACTER*8  WHICH1
      INTEGER  IS(2)
      DATA  ECCM,ECCM2  /0.002,0.00000399/
*
*
*       Define indices of KS components.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
*
*       Skip if penetration is not significant (prevents frequent calls).
      R1 = MAX(RADIUS(I1),RADIUS(I2))
      IF (KZ(27).EQ.1.AND.ABS(QPERI - 4.0*R1).LT.0.01*QPERI) GO TO 50
*
*       Skip procedure if both stars are highly evolved.
      IF (RADIUS(I1) + RADIUS(I2).LE.0.0D0) GO TO 50
*
*       Set c.m. index & reduced mass and increase event counter.
      I = N + IPAIR
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
      NDISS = NDISS + 1
      RKS = R(IPAIR)
*
*       Determine pericentre variables U & UDOT by backwards integration.
      CALL KSPERI(IPAIR)
*
*       Form semi-major axis & eccentricity.
      TD2 = 0.0
      DO 5 K = 1,4
          TD2 = TD2 + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
    5 CONTINUE
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TD2**2/(SEMI*BODY(I))
      ECC = SQRT(ECC2)
      AM0 = SEMI*(1.0D0 - ECC**2)
*
      PERI = SEMI*(1.0D0 - ECC)
      HI = H(IPAIR)
*
*       Copy stellar type.
      IS(1) = KSTAR(I1)
      IS(2) = KSTAR(I2)
*
*       Obtain kinetic energy loss due to tidal interaction (DE > 0 here).
      CALL TIDES(QPERI,BODY(I1),BODY(I2),RADIUS(I1),RADIUS(I2),IS,DE)
*
*       Include safety check on energy loss to prevent new SEMI < R.
      DH = -(DE(1) + DE(2))/ZMU
      IF (H(IPAIR) + DH.LT.-0.5*BODY(I)/R(IPAIR)) THEN
          DH = -0.5*BODY(I)/R(IPAIR) - H(IPAIR)
          DE(1) = -ZMU*DH
          DE(2) = 0.0
      END IF
*
*       Suppress the old PT procedure (DH => ECC from AM0 = const).
*     ECC2 = ECC**2 + 2.0D0*AM0*DH/BODY(I)
*     ECC2 = MAX(ECC2,ECCM2)
*
*       Adopt sequential circularization instead of standard PT.
      ECC2 = ECCM2
      ECC1 = SQRT(ECC2)
      ACIRC = AM0/(1.0 - ECC2)
*
*       Accept circularized orbit if ACIRC < 4*R1 (use maximum radius).
      IF (ACIRC.LT.4.0*R1) THEN
          SEMI1 = ACIRC
      ELSE
*       Obtain E1 from A1*(1 - E1) = 4*R1 and A1 from AM = const.
          ECC1 = 0.25*AM0/R1 - 1.0
          ECC1 = MAX(ECC1,ECCM)
          ECC1 = MAX(ECC1,0.9D0*ECC)
          SEMI1 = AM0/(1.0 - ECC1**2)
      END IF
*
*       Form the corresponding energy change.
      DH = 0.5*BODY(I)*(1.0/SEMI - 1.0/SEMI1)
*     DH = (ECC2 - ECC**2)*BODY(I)/(2.0D0*AM0)
      DE(1) = -ZMU*DH
      DE(2) = 0.0
*
*       Skip possible hyperbolic case.
      IF (H(IPAIR) + DH.GT.0.0) GO TO 50
*
*       Update total energy loss.
      ECOLL = ECOLL + (DE(1) + DE(2))
      E(10) = E(10) + (DE(1) + DE(2))
      H(IPAIR) = H(IPAIR) + DH
*
*       Print first and last energy change and update indicator.
      IF (KSTAR(I).EQ.0.OR.KSTAR(I).EQ.9) THEN
          P = DAYS*SEMI1*SQRT(SEMI1/BODY(I))
          IF (KSTAR(I).EQ.0.AND.ECC1.GT.ECCM) THEN
              WHICH1 = 'NEW CIRC'
              KSTAR(I) = 9
          ELSE
              WHICH1 = 'END CIRC'
              KSTAR(I) = 10
              NCIRC = NCIRC + 1
          END IF
          WRITE (6,8)  WHICH1, NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                 TTOT, ECC, ECC1, P, SEMI1, R1
    8     FORMAT (' ',A8,'    NAM K* T E0 EF P AF R* ',
     &                        2I6,2I4,F9.2,2F8.3,F7.1,1P,2E10.2)
      END IF
*
*       Set new pericentre.
   10 PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*
*       Specify KS velocity scaling (conserved J or chaos treatment).
      IF (KZ(27).EQ.1) THEN
          C2 = 1.0/C1
      ELSE
          C2 = SQRT((BODY(I) + H(IPAIR)*PERI1)/(BODY(I) + HI*R(IPAIR)))
      END IF
*
*       See whether circular orbit condition applies.
*     IF (ECC1.LE.ECCM) THEN
*         AM = SEMI1*(1.0D0 - ECC1**2)
*         C2 = SQRT(AM/AM0)/C1
*     END IF
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
      DO 15 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
   15 CONTINUE
*
*       Form new perturber list after significant energy loss.
      NP0 = LIST(1,I1)
      IF (ABS(SEMI1/SEMI).LT.0.5) THEN
          CALL KSLIST(IPAIR)
      END IF
*
*       Ensure chaos treatment at every pericentre by perturbed motion.
      IF (KZ(27).GT.1.AND.LIST(1,I1).EQ.0) THEN
          LIST(1,I1) = 1
          LIST(2,I1) = N
          NP0 = 1
      END IF
*
*       Re-initialize KS polynomials at pericentre for perturbed case.
      T0(I1) = TIME
      IF (NP0.GT.0) THEN
          CALL RESOLV(IPAIR,1)
          IMOD = KSLOW(IPAIR)
          CALL KSPOLY(IPAIR,IMOD)
      END IF
*
      IF (ECC.GT.0.99.AND.ABS(ECC - ECC1).GT.0.01.OR.KZ(27).GT.2) THEN
          WRITE (6,20)  NAME(I1), NAME(I2), SEMI1, ECC, ECC1, HI, QPERI
   20     FORMAT (' NEW KSTIDE    NAM AF E0 EF HI QP ',
     &                            2I5,1PE10.2,0P2F8.3,F9.1,1PE10.2)
      END IF
*
*       Set one unperturbed period for small apocentre perturbation (#27=1).
      GA = GAMMA(IPAIR)*(SEMI1*(1.0 + ECC1)/R(IPAIR))**3
      IF (KZ(27).EQ.1.AND.GA.LT.GMIN) THEN
          STEP(I1) = TWOPI*SEMI1*SQRT(SEMI1/BODY(I))
          LIST(1,I1) = 0
      END IF
*
*       Ensure T'' = 0 for pericentre test in KSINT & dissipation in UNPERT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Count any hyperbolic captures.
      IF (SEMI.LT.0.0.AND.SEMI1.GT.0.0) THEN
          NTIDE = NTIDE + 1
      END IF
*
*       Record diagnostics for new synchronous orbit and activate indicator.
      RX = MAX(RADIUS(I1),RADIUS(I2))
      IF (ECC.GT.ECCM.AND.ECC1.LT.ECCM.AND.KZ(27).EQ.1) THEN
          NSYNC = NSYNC + 1
          ESYNC = ESYNC + ZMU*H(IPAIR)
          KSTAR(I) = 10
*         WRITE (6,35)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), SEMI1,
*    &                  ECC, ECC1, HI, QPERI, RX
*  35     FORMAT (' CIRCULARIZED    NAM K* AF E0 EF HI QP R* ',
*    &                         2I6,2I3,1P,E10.2,0P,2F8.3,F9.1,1P,2E10.2)
      END IF
*
*       See whether a low-eccentricity synchronous state has been reached.
      RCOLL = 0.75*(RADIUS(I1) + RADIUS(I2))
      IF (ABS(SEMI1).LT.1.5*RCOLL.AND.ECC.LT.ECCM.AND.
     &    KSTAR(I).LT.10) THEN
          KSTAR(I) = 10
          WRITE (6,40)  ECC1, SEMI1, R(IPAIR), RCOLL, RX
   40     FORMAT (' INACTIVE PHASE    E A R RC R* ',F7.3,1P,4E10.2)
      END IF
*
*       Include warning in case of eccentricity increase (PT only).
      ECC2 = 1.0 - R(IPAIR)/SEMI1
      IF (KZ(27).EQ.1.AND.ECC2.GT.MAX(ECC,ECCM)) THEN
          WRITE (6,42)  TTOT, IPAIR, ECC2, ECC, R(IPAIR), SEMI1
   42     FORMAT (' WARNING!    E > E0    T IP E E0 R A ',
     &                                    F10.4,I4,2F8.4,1P,2E10.2)
      END IF
*
*       Reduce radius by DR/R to delay dissipation for small energy loss.
      IF (KZ(27).EQ.1.AND.DE(1)+DE(2).LT.1.0D-07*ZMU*ABS(HI)) THEN
          J1 = I1
          IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
***       IF (TEV(J1) - TIME.GT.0.01*TEV(J1)) TEV(J1) = TIME
          DR = (0.99*4.0*RADIUS(J1) - QPERI)/QPERI
          YF = MAX(ABS(DR),0.01D0)
          RADIUS(J1) = (1.0 - MIN(YF,0.1D0))*RADIUS(J1)
          DE1 = (DE(1) + DE(2))/(ZMU*ABS(H(IPAIR)))
          WRITE (6,44)  TTOT, KSTAR(J1), H(IPAIR), QPERI, DR, DE1
   44     FORMAT (' REDUCED RADIUS    T K* H QP DR/R DE/E ',
     &                                F9.2,I3,F10.1,1P,3E9.1)
      END IF
*
*       Combine the two stars inelastically in case of chaos disruption.
   45 IF (KZ(27).GT.1.AND.IDIS.GT.0) THEN
          WRITE (6,48) TTOT, IPAIR, LIST(1,I1), ECC, SEMI, QPERI,
     &                 RADIUS(I1), RADIUS(I2)
   48     FORMAT (' DISRUPTED CHAOS    T KS NP E A QP R* ',
     &                                 F9.2,2I4,F8.4,1P,4E10.2)
          CALL XVPRED(I,0)
          KSPAIR = IPAIR
          IQCOLL = 1
          CALL CMBODY(R(IPAIR),2)
      END IF
*
   50 RETURN
*
      END
