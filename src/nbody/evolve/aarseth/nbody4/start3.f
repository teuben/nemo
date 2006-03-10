      SUBROUTINE START3(ISUB)
*
*
*       Initialization & restart of triple system.
*       ------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  M
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/RCLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
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
*       Define sequence with dominant binary component as reference body.
      NAME3(1) = JCLOSE
*       Perturber index from routine IMPACT (ICOMP & JCOMP set in KSTERM).
      IF (BODY(ICOMP).LT.BODY(JCOMP)) THEN
          K = ICOMP
          ICOMP = JCOMP
          JCOMP = K
      END IF
      NAME3(2) = JCOMP
      NAME3(3) = ICOMP
*
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,3
          J = NAME3(L)
          M(L) = BODY(J)
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X3(K,L) = X(K,J)
              XDOT3(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X3(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT3(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities for triple system.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for three-body regularization.
      DO 8 L = 1,3
          DO 7 K = 1,3
              X3(K,L) = X3(K,L) - CM(K)
              XDOT3(K,L) = XDOT3(K,L) - CM(K+3)
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and include in total subsystem energy.
      CALL TRANS3(0)
      ESUB = ESUB + 0.5D0*ENERGY
*
*       Save global indices & particle attributes of subsystem.
      DO 10 L = 1,3
          J = NAME3(L)
          JLIST(L) = J
          NAME3(L) = NAME(J)
          SIZE(L) = RADIUS(J)
          IP(L) = KSTAR(J)
   10 CONTINUE
*
*       Make perturber list for potential energy correction.
      ICM = JLIST(1)
      GO TO 200
*
*       Obtain potential energy of resolved subsystem & perturbers.
   20 CALL NBPOT(3,NP,POT1)
*
*       Define subsystem indicator (ISYS = 1, 2, 3 for triple, quad, chain).
      ISYS(NSUB+1) = 1
*
*       Form ghosts and initialize c.m. motion in ICOMP (= JLIST(1) = I).
      CALL SUBSYS(3,CM)
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
      CALL NBREM(ICM,3,NPAIRS)
*
*       Set maximum integration interval equal to c.m. step.
      TMAX = STEP(ICOMP)
*
*       Copy total energy and output & capture option for routine TRIPLE.
      CM(8) = BE(3)
      KZ15 = KZ(15)
      KZ27 = KZ(27)
*
*       Assign new subsystem index and begin triple regularization.
      ISUB = NSUB
      NTRIP = NTRIP + 1
      GO TO 180
*
*       Prepare KS regularization and direct integration of third body.
  100 IMIN = 1
      IF (R2.LT.R1) IMIN = 2
*
*       Specify global names of the KS candidates and least dominant body.
      NAM1 = NAME3(3)
      NAM2 = NAME3(IMIN)
      NAM3 = NAME3(3-IMIN)
*
*       Identify current global indices by searching all single particles.
      DO 102 J = IFIRST,N
          IF (NAME(J).EQ.NAM1) ICOMP = J
          IF (NAME(J).EQ.NAM2) JCOMP = J
          IF (NAME(J).EQ.NAM3) I = J
  102 CONTINUE
*
*       Identify global index for c.m. body.
      ICM = I
      IF (BODY(ICOMP).GT.0.0D0) ICM = ICOMP
      IF (BODY(JCOMP).GT.0.0D0) ICM = JCOMP
*
*       Distinguish between block-step scheme and standard integration.
      IF (TBLOCK.GT.0.0) THEN
*       Define discrete time for new polynomials (subject to TIME <= TBLOCK).
          DT = 0.1*MIN(TIME3,STEP(ICM))
          CALL STEPK(DT,DTN)
          TIME2 = T0S(ISUB) + TIME3 - TPREV
          TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
          TIME = MIN(TBLOCK,TIME)
          IF (DMOD(TIME,DTK(40)).NE.0.0D0) TIME = TBLOCK
      ELSE
*       Update the global time and delay output times to avoid troubles.
          TIME = T0S(ISUB) + TIME3
          TADJ = MAX(TADJ,TIME)
          TPRINT = MAX(TPRINT,TIME)
      END IF
*
*       Predict current coordinates & velocities to F3DOT before termination.
      CALL XVPRED(ICM,-1)
*
*       Copy c.m. coordinates & velocities.
      DO 105 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
  105 CONTINUE
*
*       Re-determine the perturber list.
      GO TO 200
*
*       Obtain potential energy of the c.m. subsystem & JPERT(NP).
  110 JLIST(1) = ICM
      CALL NBPOT(1,NP,POT1)
*
*       Save global indices of three-body system.
      JLIST(1) = I
      JLIST(2) = ICOMP
      JLIST(3) = JCOMP
*
*       Set configuration pointers for escaper & KS candidates.
      JLIST(4) = 3 - IMIN
      JLIST(5) = 3
      JLIST(6) = IMIN
*
*       Place new coordinates in the original locations.
      DO 120 L = 1,3
          J = JLIST(L)
*       Compare global name & subsystem name to restore the mass.
          DO 112 K = 1,3
              IF (NAME(J).EQ.NAMES(K,ISUB)) THEN
                  BODY(J) = BODYS(K,ISUB)
              END IF
  112     CONTINUE
          LL = JLIST(L+3)
          DO 115 K = 1,3
              X(K,J) = X3(K,LL) + CM(K)
  115     CONTINUE
  120 CONTINUE
*
*       Obtain potential energy of subsystem & perturbers at the end.
      CALL NBPOT(3,NP,POT2)
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
      DO 130 L = 1,3
          J = JLIST(L)
          LL = JLIST(L+3)
          DO 125 K = 1,3
              XDOT(K,J) = XDOT3(K,LL) + CM(K+3)
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
      ESUB = ESUB - 0.5D0*ENERGY - ECOLL3
*
*       Replace ICM in perturber lists by all subsystem members.
      CALL NBREST(ICM,3,NPAIRS)
*
*       Check for star collision (only needs coordinates & velocities).
      IF (ITERM.LT.0) THEN
          JLIST(1) = MIN(ICOMP,JCOMP)
          JLIST(2) = MAX(ICOMP,JCOMP)
          JLIST(3) = I
          JLIST(4) = 0
          GO TO 170
      END IF
*
*       Initialize force polynomial and time-step for body #I.
      CALL FPOLY1(I,I,0)
*
*       See whether body #I should be added to NLIST.
      IF (T0(I) + STEP(I).LT.TLIST) THEN
          CALL NLMOD(I,1)
      END IF
*
*       Perform KS regularization of dominant components (ICOMP < JCOMP).
      I = ICOMP
      ICOMP = MIN(ICOMP,JCOMP)
      JCOMP = MAX(I,JCOMP)
      CALL KSREG
*
*       Check minimum two-body distance.
      DMIN3 = MIN(DMIN3,RCOLL)
*
*       Update net binary energy change.
      SBCOLL = SBCOLL + CM(9)
*
*       Set IPHASE = -1 since routine INTGRT exits on IPHASE > 0.
      IPHASE = -1
*
*       Update number of DIFSY calls, tidal dissipations & collision energy.
  170 NSTEPT = NSTEPT + NSTEP3
      NDISS = NDISS + NDISS3
      ECOLL = ECOLL + ECOLL3
      E(10) = E(10) + ECOLL3
*
*       Check for subsystem at last COMMON dump (no restart with NSUB > 0).
      IF (NSUB.EQ.0.AND.KZ(2).GE.1) THEN
          IF (TIME - TDUMP.LT.TIME3) THEN
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
              IF (J.EQ.ICM) GO TO 210
              IF (J.EQ.ICOMP.OR.J.EQ.JCOMP) GO TO 210
              NP = NP + 1
              JPERT(NP) = J
          END IF
  210 CONTINUE
*
      IF (ISUB.EQ.0) GO TO 20
      GO TO 110
*
      END
