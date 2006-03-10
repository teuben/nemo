      SUBROUTINE SETUP
*
*
*       Generation of initial coordinates & velocities.
*       -----------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(8)
      REAL*4  RAN2
*
*
*       Choose between uniform density and Plummer model.
      KDUM = IDUM1
      IF (KZ(5).GT.0) GO TO 20
*
*       Set up a uniform spherical system.
      DO 10 I = 1,N
    1     A(1) = 0.0D0
          DO 2 K = 1,3
              A(K+1) = 2.0*RAN2(KDUM) - 1.0
              A(1) = A(1) + A(K+1)**2
    2     CONTINUE
          IF (A(1).GT.1.0) GO TO 1
*
    4     A(5) = 0.0D0
          DO 5 K = 1,3
              A(K+5) = 2.0*RAN2(KDUM) - 1.0
              A(5) = A(5) + A(K+5)**2
    5     CONTINUE
          IF (A(5).GT.1.0) GO TO 4
*
          DO 8 K = 1,3
*             X(K,I) = A(1)*A(K+1)
*       Density proportional to 1/R**2.
              X(K,I) = A(K+1)
*       Constant density.
              XDOT(K,I) = A(K+5)
*       Isotropic velocities (magnitude randomized; no radial dependence).
    8     CONTINUE
   10 CONTINUE
*
      GO TO 90
*
*       Initialize centre of mass terms.
   20 DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
*       Generate initial conditions from Plummer model (A & A 37, 183).
      DO 40 I = 1,N
   30     A(1) = RAN2(KDUM)
          IF (A(1).LT.1.0D-10) GO TO 30
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 30
*
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
   32     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 32
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*
*       Accumulate c.m. terms.
          DO 35 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   35     CONTINUE
   40 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 50 I = 1,N
          DO 45 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
   45     CONTINUE
*          IF(KZ(22).EQ.1) WRITE(10,*)(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
   50 CONTINUE
      IF(KZ(22).EQ.1) CALL FLUSH(10)
*
*       Check initial conditions for two orbiting Plummer spheres.
      IF (KZ(5).LE.1) GO TO 90
*
*       Save membership of first system for colour plot (N2 = NZERO - N1).
      N1 = N
      IF (KZ(5).EQ.2) THEN
      READ (5,*)  APO, ECC, N2, SCALE
      N2 = MIN(N,N2)
      SEMI = APO/(1.0 + ECC)
      SEMI = MIN(SEMI,50.0D0)
      SEMI = MAX(SEMI,2.0D0)
      ECC = MIN(ECC,0.999D0)
      ECC = MAX(ECC,0.0D0)
      ZM2 = 0.0
      KSKIP = N1/N2
      DO 52 I = 1,N2
          J = KSKIP*I
          ZM2 = ZM2 + BODY(J)
   52 CONTINUE
      FAC1 = ZM2/(ZMASS + ZM2)
      FAC2 = ZMASS/(ZMASS + ZM2)
*       Restrict volume ratio to 125 (i.e. unreasonable density contrast).
      IF (SCALE.LE.0.2D0) SCALE = 0.2
*       Increase total mass (save in ZMTOT for possible use in XTRNL0).
      ZMTOT = ZMTOT + ZMBAR*ZM2
      ZMASS = ZMASS + ZM2
*       Set apocentre velocity for new combined mass.
      VAP = SQRT(ZMASS/SEMI)*SQRT((1.0 - ECC)/(1.0 + ECC))
      DO 55 I = 1,N
          IF (I.LE.N2) THEN
*       Copy members from first system by uniform skipping (N2 <= N1).
          J = KSKIP*I
          BODY(I+N) = BODY(J)
          X(1,I+N) = SCALE*X(1,J) + FAC2*APO
          X(2,I+N) = SCALE*X(2,J)
          X(3,I+N) = SCALE*X(3,J)
          XDOT(1,I+N) = XDOT(1,J)/SQRT(SCALE)
          XDOT(2,I+N) = XDOT(2,J)/SQRT(SCALE) + FAC2*VAP
          XDOT(3,I+N) = XDOT(3,J)/SQRT(SCALE)
          END IF
          X(1,I) = X(1,I) - FAC1*APO
          XDOT(2,I) = XDOT(2,I) - FAC1*VAP
   55 CONTINUE
      ELSE IF (KZ(5).EQ.3) THEN
*       Prepare case of accretion disk with massive perturber.
          READ (5,*)  APO, ECC, DMIN, SCALE
          RIN = 0.5
          ROUT = 1.0
          ZMASS = 1.0
          BODY(1) = ZMASS
          DO 99 K = 1,3
              X(K,1) = 0.0
              XDOT(K,1) = 0.0
   99     CONTINUE
*       Generate a thin disk population in circular orbits.
          DO 100 I = 2,N
              BODY(I) = 1.0D-03/FLOAT(N)
              SEMI = RIN + (ROUT - RIN)*FLOAT(I)/FLOAT(N)
              VCIRC = SQRT((BODY(1) + BODY(I))/SEMI)
              PHASE = TWOPI*RAN2(KDUM)
              X(1,I) = SEMI*COS(PHASE)
              X(2,I) = SEMI*SIN(PHASE)
              X(3,I) = 0.01*(2.0*RAN2(KDUM) - 1.0)
              XDOT(1,I) = -VCIRC*SIN(PHASE)
              XDOT(2,I) = VCIRC*COS(PHASE)
              XDOT(3,I) = 0.01*VCIRC*(2.0*RAN2(KDUM) - 1.0)
  100     CONTINUE
*       Define membership of perturber and ensure no external tide.
          N2 = 1
          KZ(14) = 0
*       Redefine solar mass unit and astronomical length scale in AU.
          ZMBAR = 1.0
          RBAR = 1.0/2.05D+05
          BODY(N+1) = SCALE*BODY(1)
          ZMTOT = ZMASS + BODY(N+1)
*       Set appropriate mass ratios for transforming to new c.m. frame.
          FAC1 = BODY(N+1)/(ZMASS + BODY(N+1))
          FAC2 = ZMASS/(ZMASS + BODY(N+1))
          ZMASS = ZMASS + BODY(N+1)
*       Form orbital elements for massive perturber (avoid ECC = 1).
          IF (ABS(ECC - 1.0).GT.1.0D-05) THEN
              SEMI = DMIN/(1.0 - ECC)
          ELSE
              SEMI = -1.0D+05
          END IF
          VM2 = ZMASS*(2.0/DMIN - 1.0/SEMI)
          VAP2 = ZMASS*(2.0/APO - 1.0/SEMI)
*       Determine initial y-velocity from angular momentum conservation.
          VY = SQRT(VM2)*DMIN/APO
          VX = SQRT(VAP2 - VY**2)
*       Place perturber on the Y-axis with appropriate velocities.
          X(1,N+1) = APO*FAC2
          X(2,N+1) = 0.0
          X(3,N+1) = 0.0
          XDOT(1,N+1) = -VX*FAC2
          XDOT(2,N+1) = VY*FAC2
          XDOT(3,N+1) = 0.0
*       Displace the disk members and include negative y-velocity.
          DO 120 I = 1,N
              X(1,I) = X(1,I) - FAC1*APO
              XDOT(1,I) = XDOT(1,I) + FAC1*VX
              XDOT(2,I) = XDOT(2,I) - FAC1*VY
  120     CONTINUE
      ELSE IF (KZ(5).EQ.4) THEN
*       Include two massive bodies (ECC > 1: NAME = 1 & 2 free floating).
          N2 = 0
          READ (5,*)  SEMI, ECC, ZM1, ZM2
          WRITE (6,125)  SEMI, ECC, ZM1, ZM2
  125     FORMAT (/,12X,'MASSIVE BODIES    A =',1P,E9.1,
     &            '  E =',0P,F6.2,'  M1/<M> =',F6.2,'  M2/<M> =',F6.2)
          BODY(1) = ZM1
          BODY(2) = ZM2
          IF (ECC.LT.1.0) THEN
*       Set apocentre velocity for new combined mass (using NAME = 1 & 2).
              VAP = SQRT((ZM1 + ZM2)/SEMI)*SQRT((1.0 - ECC)/(1.0 + ECC))
              FAC1 = ZM2/(ZM1 + ZM2)
              FAC2 = ZM1/(ZM1 + ZM2)
              DO 130 K = 1,3
                  X(K,1) = 0.0
                  X(K,2) = 0.0
                  XDOT(K,2) = 0.0
  130         CONTINUE
*       Initialize binary with c.m. at rest (elements change in SCALE).
              X(1,1) = -FAC1*SEMI*(1.0 + ECC)
              X(1,2) = FAC2*SEMI*(1.0 + ECC)
              XDOT(2,1) = -FAC1*VAP
              XDOT(2,2) = FAC2*VAP
          END IF
      ELSE
          GO TO 90
      END IF
*
*       Specify new membership.
      N = N + N2
      NZERO = N
      NTOT = N
      IF (N.GE.NMAX-10) THEN
          WRITE (6,56)  N, NMAX
   56     FORMAT (' DANGER!    LIMIT EXCEEDED   N =',I6,'  NMAX =',I6)
          STOP
      END IF
*
      IF (KZ(5).EQ.2) THEN
          WRITE (6,58)  SEMI, ECC, N1, N2, SCALE
   58     FORMAT (/,12X,'PLUMMER BINARY    A =',F6.2,'  E =',F6.2,
     &                  '  N1 =',I6,'  N2 =',I6,'  SCALE =',F6.2)
      ELSE IF (KZ(5).EQ.3) THEN
          WRITE (6,59)  APO, ECC, DMIN, SCALE
   59     FORMAT (/,12X,'MASSIVE PERTURBER    APO =',F6.2,'  E =',F6.2,
     &                  '  DMIN =',F6.2,'  MP/M1 =',F6.2)
      END IF
*
*       Re-initialize centre of mass terms.
      DO 60 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   60 CONTINUE
      ZMASS = 0.0
      DO 70 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 65 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   65     CONTINUE
   70 CONTINUE
      DO 80 I = 1,N
          DO 75 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   75     CONTINUE
   80 CONTINUE
*
*       Save random number sequence in COMMON for future use.
   90 IDUM1 = KDUM
*
      RETURN
*
      END
