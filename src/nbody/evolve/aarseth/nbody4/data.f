      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common4.h'
      REAL*4  RAN2
*
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 1 K = 1,IDUM1
          RN1 = RAN2(KDUM)
    1 CONTINUE
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
*       Set provisional total mass (rescaled in routine SCALE).
      ZMASS = FLOAT(N)
*
*       Read mass function parameters, # primordials, Z-abundance, epoch
*       and interval for calling routine HRPLOT.
      READ (5,*)  ALPHA, BODY1, BODYN, NBIN0, ZMET, EPOCH0, DTPLOT
      IF (KZ(25).GT.0) DTPLOT = MAX(DTPLOT,DELTAT)
*
*       Check option for reading initial conditions from input file.
      IF (KZ(22).GE.2) THEN
          ZMASS = 0.0D0
          DO 5 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
              ZMASS = ZMASS + BODY(I)
    5     CONTINUE
*       Allow own IMF and standard scaling.
          IF (KZ(22).EQ.2.AND.BODY1.EQ.BODYN) GO TO 30
*       Specify total mass as N times mean mass.
          ZMTOT = FLOAT(N)*ZMBAR
          IF (KZ(22).GE.3) GO TO 40
      END IF
*
*       Include the case of equal masses (ALPHA = 1 or BODY1 = BODYN).
      IF (ALPHA.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
*       Set provisional total mass (rescaled in routine SCALE).
          ZMTOT = FLOAT(N)*ZMBAR
          GO TO 40
      END IF
*
      CALL IMF(BODY1,BODYN)
*
*       Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR
   35 CONTINUE
      ZMTOT = ZMASS
      ZMASS = FLOAT(N)
*
*       Set up initial coordinates & velocities (uniform or Plummer model).
   40 IF (KZ(22).LE.1) THEN
          CALL SETUP
      END IF
*
   50 RETURN
*
      END
