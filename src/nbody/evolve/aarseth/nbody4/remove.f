      SUBROUTINE REMOVE(I,KCASE)
*
*
*       Particle removal.
*       -----------------
*
      INCLUDE 'common4.h'
      REAL*8  A(6)
*
*
*       Remove escaper, KS pair, components or c.m. (KCASE = 1, 2, 3).
      IF (KCASE.EQ.2) GO TO 20
      IF (KCASE.EQ.3) GO TO 10
*
*       Correct force & first derivative of all bodies (only for escape).
      DO 5 J = IFIRST,NTOT
          IF (J.EQ.I) GO TO 5
          RIJ2 = 0.0D0
          A7 = 0.0D0
          DO 1 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              A7 = A7 + A(K)*A(K+3)
    1     CONTINUE
          A8 = BODY(I)/(RIJ2*SQRT(RIJ2))
          DO 2 K = 1,3
              A(K+3) = (A(K+3) - 3.0*A7*A(K)/RIJ2)*A8
              F(K,J) = F(K,J) - 0.5*A(K)*A8
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
    2     CONTINUE
    5 CONTINUE
*
*       Move up all COMMON variables (escaper or old c.m. & KS comps).
   10 IF (I.GT.NTOT) GO TO 50
*
      DO 15 J = I,NTOT
          J1 = J + 1
          DO 12 K = 1,3
              X(K,J) = X(K,J1)
              X0(K,J) = X0(K,J1)
              X0DOT(K,J) = X0DOT(K,J1)
              XDOT(K,J) = XDOT(K,J1)
              F(K,J) = F(K,J1)
              FDOT(K,J) = FDOT(K,J1)
              D2(K,J) = D2(K,J1)
              D3(K,J) = D3(K,J1)
   12     CONTINUE
*
          BODY(J) = BODY(J1)
          RADIUS(J) = RADIUS(J1)
          ZLMSTY(J) = ZLMSTY(J1)
          SPIN(J) = SPIN(J1)
          TEV(J) = TEV(J1)
          TEV0(J) = TEV0(J1)
          BODY0(J) = BODY0(J1)
          EPOCH(J) = EPOCH(J1)
          KSTAR(J) = KSTAR(J1)
          NAME(J) = NAME(J1)
          STEP(J) = STEP(J1)
          T0(J) = T0(J1)
          PHI(J) = PHI(J1)
   15 CONTINUE
*
      GO TO 50
*
*       Move up all tables of KS pairs below IPAIR = I.
   20 DO 30 JPAIR = I,NPAIRS
          JP1 = JPAIR + 1
          DO 25 K = 1,4
              U(K,JPAIR) = U(K,JP1)
              U0(K,JPAIR) = U0(K,JP1)
              UDOT(K,JPAIR) = UDOT(K,JP1)
              FU(K,JPAIR) = FU(K,JP1)
              FUDOT(K,JPAIR) = FUDOT(K,JP1)
              FUDOT2(K,JPAIR) = FUDOT2(K,JP1)
              FUDOT3(K,JPAIR) = FUDOT3(K,JP1)
              SF(K,JPAIR) = SF(K,JP1)
              FP0(K,JPAIR) = FP0(K,JP1)
              FD0(K,JPAIR) = FD0(K,JP1)
   25     CONTINUE
*
          R(JPAIR) = R(JP1)
          R0(JPAIR) = R0(JP1)
          DTAU(JPAIR) = DTAU(JP1)
          TDOT2(JPAIR) = TDOT2(JP1)
          TDOT3(JPAIR) = TDOT3(JP1)
          GAMMA(JPAIR) = GAMMA(JP1)
          H(JPAIR) = H(JP1)
          HDOT(JPAIR) = HDOT(JP1)
          HDOT2(JPAIR) = HDOT2(JP1)
          HDOT3(JPAIR) = HDOT3(JP1)
          HDOT4(JPAIR) = HDOT4(JP1)
          KSLOW(JPAIR) = KSLOW(JP1)
          SF(5,JPAIR) = SF(5,JP1)
          SF(6,JPAIR) = SF(6,JP1)
          SF(7,JPAIR) = SF(7,JP1)
          H0(JPAIR) = H0(JP1)
*
*       Transfer two perturber lists for each pair.
          DO 29 KCOMP = 1,2
              NNB = LIST(1,2*JPAIR+KCOMP) + 1
*       Include flag of 2nd component.
              IF (NNB.EQ.1) NNB = 2
              DO 28 L = 1,NNB
                  LIST(L,2*JPAIR-2+KCOMP) = LIST(L,2*JPAIR+KCOMP)
   28         CONTINUE
   29     CONTINUE
   30 CONTINUE
*
   50 RETURN
*
      END
