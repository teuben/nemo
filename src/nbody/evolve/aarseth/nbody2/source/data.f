      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common2.h'
      REAL*4  RAN2,A(8)
*
*
*       Check option for reading initial conditions from input file.
      IF (KZ(4).EQ.2) THEN
*       Note that initial conditions are assumed REAL*4 here.
          DO 2 I = 1,N
              READ (4)  BODY(I), (A(K),K=1,3), (XDOT(K,I),K=1,3)
              DO 2 K = 1,3
                  X(K,I) = A(K)
    1         CONTINUE
    2     CONTINUE
          GO TO 90
      END IF
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      A(1) = RAN2(KDUM)
*       Skip the first random numbers (JLIST(1) specified at input).
      NRAND = JLIST(1)
      DO 5 K = 1,NRAND
          A(1) = RAN2(KDUM)
    5 CONTINUE
*
*       Read mass parameters.
      READ (5,*)  ALPHAS, BODY1, BODYN
*
*       Include the case of equal masses (ALPHAS = 1 or BODY1 = BODYN).
      IF (ALPHAS.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
          ZMASS = FLOAT(N)
          GO TO 25
      END IF
*
      WRITE (6,15)  ALPHAS, BODY1, BODYN
   15 FORMAT (/,12X,'STANDARD MASS FUNCTION:','   ALPHAS =',F5.2,
     &                                '  BODY1 =',F5.1,'  BODYN =',F5.2)
*
*       Generate a power-law mass function with exponent ALPHAS.
      ALPHA1 = ALPHAS - 1.0
      FM1 = 1.0/BODY1**ALPHA1
      FMN = (FM1 - 1.0/BODYN**ALPHA1)/(FLOAT(N) - 1.0)
      ZMASS = 0.0
      CONST = 1.0/ALPHA1
*
      DO 18 I = 1,N
          FMI = FM1 - FLOAT(I - 1)*FMN
          BODY(I) = 1.0/FMI**CONST
          ZMASS = ZMASS + BODY(I)
   18 CONTINUE
*
*       First scale the masses to <M> = 1.
      ZMBAR1 = ZMASS/FLOAT(N)
      DO 20 I = 1,N
          BODY(I) = BODY(I)/ZMBAR1
   20 CONTINUE
      ZMASS = FLOAT(N)
*
   25 IF (KZ(5).GT.0) GO TO 80
*
*       Set up a uniform spherical system.
      DO 40 I = 1,N
   32     A(1) = 0.0
          DO 33 K = 1,3
              A(K+1) = 2.0*RAN2(KDUM) - 1.0
              A(1) = A(1) + A(K+1)**2
   33     CONTINUE
          IF (A(1).GT.1.0) GO TO 32
   34     A(5) = 0.0
          DO 35 K = 1,3
              A(K+5) = 2.0*RAN2(KDUM) - 1.0
              A(5) = A(5) + A(K+5)**2
   35     CONTINUE
          IF (A(5).GT.1.0) GO TO 34
          DO 36 K = 1,3
*             X(K,I) = A(1)*A(K+1)
*       Density proportional to 1/R**2.
              X(K,I) = A(K+1)
*       Constant density.
              XDOT(K,I) = A(K+5)
*       Isotropic velocities (magnitude randomized; no radial dependence).
   36     CONTINUE
   40 CONTINUE
*
      GO TO 90
*
   80 DO 81 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   81 CONTINUE
*
*       Generate initial conditions from Plummer model (A & A 37, 183).
      TWOPI = 8.0*ATAN(1.0D0)
      DO 85 I = 1,N
   82     A(1) = RAN2(KDUM)
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 82
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
   83     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 83
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*
*       Accumulate the centre of mass terms.
          DO 84 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   84     CONTINUE
   85 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 88 I = 1,N
          DO 86 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
   86     CONTINUE
   88 CONTINUE
*
   90 RETURN
*
      END
