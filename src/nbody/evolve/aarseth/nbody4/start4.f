      SUBROUTINE START4(ISUB)
*
*
*       Initialization & restart of four-body system.
*       ---------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  M,MIJ
      COMMON/CREG/  M(4),X4(3,4),XDOT4(3,4),P(12),Q(12),TIME4,ENERGY,
     &              EPSR2,XR(9),W(9),RR(6),TA(6),MIJ(6),CM(10),RMAX4,
     &              TMAX,DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/CONFIG/  R2(4,4),I1,I2,I3,I4
      COMMON/RCLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Decide between new run, termination or collision (= 0, > 0, < 0).
      IF (ISUB.NE.0) THEN
          ITERM = ISUB
          ISUB = IABS(ISUB)
          GO TO 100
      END IF
*
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,4
          J = 2*NPAIRS + L
          JLIST(L) = J
*       Copy four-body system from the first single particle locations.
          M(L) = BODY(J)
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X4(K,L) = X(K,J)
              XDOT4(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities of four-body system.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for chain regularization.
      DO 8 L = 1,4
          DO 7 K = 1,3
              X4(K,L) = X4(K,L) - CM(K)
              XDOT4(K,L) = XDOT4(K,L) - CM(K+3)
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and include in total subsystem energy.
      CALL NEWSYS(X4,XDOT4,M,4,ENERGY,GAM)
      ESUB = ESUB + ENERGY
*
*       Save global indices & particle attributes of subsystem.
      DO 10 L = 1,4
          J = JLIST(L)
          NAME4(L) = NAME(J)
          SIZE(L) = RADIUS(J)
          IP(L) = KSTAR(J)
   10 CONTINUE
*
*       Make perturber list for potential energy correction.
      ICM = JLIST(1)
      GO TO 200
*
*       Obtain potential enegy of resolved subsystem & perturbers.
   20 CALL NBPOT(4,NP,POT1)
*
*       Define subsystem indicator (ISYS = 1, 2, 3 for triple, quad, chain).
      ISYS(NSUB+1) = 2
*
*       Form ghosts and initialize c.m. motion in ICOMP (= JLIST(1)).
      CALL SUBSYS(4,CM)
*
*       Include interaction of subsystem c.m. & perturbers for net effect.
      CALL NBPOT(1,NP,POT2)
*
*       Form square of c.m. velocity correction due to differential force.
      VI2 = X0DOT(1,ICOMP)**2 + X0DOT(2,ICOMP)**2 + X0DOT(3,ICOMP)**2
      CORR = 1.0 + 2.0*(POT2 - POT1)/(CM(7)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Modify c.m. velocity by net tidal energy correction.
      DO 30 K = 1,3
          X0DOT(K,ICOMP) = SQRT(CORR)*X0DOT(K,ICOMP)
   30 CONTINUE
*
*       Remove ghosts from all perturber lists.
      CALL NBREM(ICM,4,NPAIRS)
*
*       Set maximum integration interval equal to c.m. step.
      TMAX = STEP(ICOMP)
*
*       Copy total energy and output & capture options for routine QUAD.
      CM(8) = BE(3)
      KZ15 = KZ(15)
      KZ27 = KZ(27)
*
*       Assign new subsystem index and begin four-body regularization.
      ISUB = NSUB
      NQUAD = NQUAD + 1
      GO TO 180
*
*       Prepare KS regularization & direct integration of two bodies.
  100 JLIST(5) = NAME4(I1)
      JLIST(6) = NAME4(I2)
      JLIST(7) = NAME4(I3)
      JLIST(8) = NAME4(I4)
*
*       Identify current global index by searching all single particles.
      DO 102 J = IFIRST,N
          DO 101 L = 1,4
              IF (NAME(J).EQ.JLIST(L+4)) THEN
                  JLIST(L) = J
              END IF
  101     CONTINUE
  102 CONTINUE
*
*       Ensure ICOMP < JCOMP for KS regularization.
      ICOMP = MIN(JLIST(1),JLIST(2))
      JCOMP = MAX(JLIST(1),JLIST(2))
*
*       Identify global index of c.m. body.
      DO 104 L = 1,4
          J = JLIST(L)
          IF (BODY(J).GT.0.0D0) ICM = J
  104 CONTINUE
*
*       Distinguish between block-step scheme and standard integration.
      IF (TBLOCK.GT.0.0) THEN
*       Define discrete time for new polynomials (subject to TIME <= TBLOCK).
          DT = 0.1*MIN(TIME4,STEP(ICM))
          CALL STEPK(DT,DTN)
          TIME2 = T0S(ISUB) + TIME4 - TPREV
          TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
          TIME = MIN(TBLOCK,TIME)
          IF (DMOD(TIME,DTK(40)).NE.0.0D0) TIME = TBLOCK
      ELSE
*       Update the global time and delay output times to avoid troubles.
          TIME = T0S(ISUB) + TIME4
          TADJ = MAX(TADJ,TIME)
          TPRINT = MAX(TPRINT,TIME)
      END IF
*
*       Predict current coordinates & velocities to F3DOT before termination.
      CALL XVPRED(ICM,-1)
*
*       Coopy c.m. coordinates & velocities.
      DO 105 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
  105 CONTINUE
*
*       Re-determine the perturber list.
      GO TO 200
*
*       Obtain potential energy of the c.m. subsystem & JPERT(NP).
  110 I = JLIST(1)
      JLIST(1) = ICM
      CALL NBPOT(1,NP,POT1)
*
*       Set configuration pointers for KS candidates & distant bodies.
      JLIST(5) = I1
      JLIST(6) = I2
      JLIST(7) = I3
      JLIST(8) = I4
*
*       Place new coordinates in the original locations.
      JLIST(1) = I
      DO 120 L = 1,4
          J = JLIST(L)
*       Compare global name & subsystem name to restore the mass.
          DO 112 K = 1,4
              IF (NAME(J).EQ.NAMES(K,ISUB)) THEN
                  BODY(J) = BODYS(K,ISUB)
              END IF
  112     CONTINUE
          LL = JLIST(L+4)
          DO 115 K = 1,3
              X(K,J) = X4(K,LL) + CM(K)
  115     CONTINUE
  120 CONTINUE
*
*       Obtain potential energy of subsystem & perturbers at the end.
      CALL NBPOT(4,NP,POT2)
*
*       Form square of c.m. velocity correction due to differential force.
      VI2 = CM(4)**2 + CM(5)**2 + CM(6)**2
      CORR = 1.0 + 2.0*(POT2 - POT1)/(CM(7)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Modify c.m. velocity by net tidal energy correction.
      DO 122 K = 1,3
          CM(K+3) = SQRT(CORR)*CM(K+3)
  122 CONTINUE
*
*       Transform to global velocities using corrected c.m. velocity.
      DO 130 L = 1,4
          J = JLIST(L)
          LL = JLIST(L+4)
          DO 125 K = 1,3
              XDOT(K,J) = XDOT4(K,LL) + CM(K+3)
              X0DOT(K,J) = XDOT(K,J)
  125     CONTINUE
  130 CONTINUE
*
*       Predict coordinates & velocities of perturbers to order FDOT.
      DO 140 L = 1,NP
          J = JPERT(L)
          CALL XVPRED(J,0)
*       Reduce time-step and check NLIST membership.
*         STEP(J) = MAX(0.5D0*STEP(J),TIME - T0(J))
*         IF (T0(J) + STEP(J).LT.TLIST) THEN
*             CALL NLMOD(J,1)
*         END IF
  140 CONTINUE
*
*       Update subsystem COMMON variables unless last or only case.
      IF (ISUB.LT.NSUB) THEN
          DO 150 L = ISUB,NSUB
              DO 145 K = 1,4
                  BODYS(K,L) = BODYS(K,L+1)
                  NAMES(K,L) = NAMES(K,L+1)
  145         CONTINUE
              T0S(L) = T0S(L+1)
              TS(L) = TS(L+1)
              STEPS(L) = STEPS(L+1)
              RMAXS(L) = RMAXS(L+1)
              ISYS(L) = ISYS(L+1)
  150     CONTINUE
      END IF
*
*       Reduce subsystem counter and subtract internal binding energy.
      NSUB = NSUB - 1
      ESUB = ESUB - ENERGY - ECOLL3
*
*       Replace ICM in perturber lists by all subsystem members.
      CALL NBREST(ICM,4,NPAIRS)
*
*       Check for star collision (only needs coordinates & velocities).
      IF (ITERM.LT.0) THEN
          JLIST(1) = ICOMP
          JLIST(2) = JCOMP
*
*       See whether relabelling is required.
          IF (R2(I1,I4).LT.R2(I1,I3).OR.R2(I3,I4).LT.R2(I1,I3)) THEN
              IF (R2(I1,I4).LT.R2(I3,I4)) THEN
*       Switch body #I3 & I4 to give new dominant pair I1 & I3.
                  I = JLIST(4)
                  JLIST(4) = JLIST(3)
                  JLIST(3) = I
              ELSE
*       Set JLIST(5) < 0 to denote that body #I3 & I4 will be new KS pair.
                  JLIST(5) = -1
              END IF
          END IF
          GO TO 170
      END IF
*
*       Save global indices of least dominant bodies.
      I3 = JLIST(3)
      I4 = JLIST(4)
*
*       Initialize force polynomials and time-steps for body #I3 & #I4.
      CALL FPOLY1(I3,I3,0)
      CALL FPOLY1(I4,I4,0)
*
*       See whether body #I3 or #I4 should be added to NLIST.
      I = I3
  160 IF (T0(I) + STEP(I).LT.TLIST) THEN
          CALL NLMOD(I,1)
      END IF
      IF (I.EQ.I3) THEN
          I = I4
          GO TO 160
      END IF
*
*       Perform KS regularization of dominant components (ICOMP < JCOMP).
      CALL KSREG
*
*       Check minimum two-body distance.
      DMIN4 = MIN(DMIN4,RCOLL)
*
*       Update net binary energy change.
      BBCOLL = BBCOLL + CM(9)
*
*       Set IPHASE = -1 since routine INTGRT exits on IPHASE > 0.
      IPHASE = -1
*
*       Update number of DIFSY calls, tidal dissipations & collision energy.
  170 NSTEPQ = NSTEPQ + NSTEP4
      NDISS = NDISS + NDISS4
      ECOLL = ECOLL + ECOLL3
      E(10) = E(10) + ECOLL3
*
*       Check for subsystem at last COMMON dump (no restart with NSUB > 0).
      IF (NSUB.EQ.0.AND.KZ(2).GE.1) THEN
          IF (TIME - TDUMP.LT.TIME4) THEN
              TDUMP = TIME
              CALL MYDUMP(1,2)
          END IF
      END IF
*
  180 RETURN
*
*       Form the current perturber list.
  200 RP2 = CMSEP2*RMIN**2
      NP = 0
*
*       Loop over all single particles & c.m. but skip subsystem members.
      DO 210 J = IFIRST,NTOT
          RIJ2 = (X(1,J) - CM(1))**2 + (X(2,J) - CM(2))**2 +
     &                                 (X(3,J) - CM(3))**2
          IF (RIJ2.LT.RP2) THEN
              DO 205 K = 1,4
                  IF (J.EQ.JLIST(K)) GO TO 210
  205         CONTINUE
              NP = NP + 1
              JPERT(NP) = J
          END IF
  210 CONTINUE
*
      IF (ISUB.EQ.0) GO TO 20
      GO TO 110
*
      END
