      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'commonp.h'
      REAL*4  RAN2
*
*
*       Read initial conditions.
      DO 1 I = 1,N
          READ (5,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
    1 CONTINUE
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 10 K = 1,IDUM1
          RN1 = RAN2(KDUM)
   10 CONTINUE
*
*       Initialize centre of mass terms.
      DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
*     DO 40 I = 1,NMASS
*         ZMASS = ZMASS + BODY(I)
*         DO 35 K = 1,3
*             XDOT(K,I) = VSCALE*XDOT(K,I)
*             CMR(K) = CMR(K) + BODY(I)*X(K,I)
*             CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
*  35     CONTINUE
*  40 CONTINUE
*
*       Define nominal crossing time (cf. initial step and zero velocity).
      TCR = TWOPI
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
      RETURN
*
      END
