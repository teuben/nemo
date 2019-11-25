      SUBROUTINE START
*
*
*       Polynomial initialization.
*       --------------------------
*
      INCLUDE 'common1.h'
      EXTERNAL SCALE
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
*
*       Check for optional creation of two separate subsystems.
      IF (KZ(8).GT.0) THEN
          CALL SUBSYS
      END IF
*
*       Set sequential name & total mass.
      ZMASS = 0.0
      DO 10 I = 1,N
          NAME(I) = I
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
*       Obtain force & first derivative.
      CALL FPOLY1(1,N)
*
*       Obtain second & third force derivatives.
      CALL FPOLY2(1,N)
*
*       Set new time-steps and initialize divided differences.
      CALL STEPS(1,N)
*
*       Initialize the time-step list used to find next body.
      DTLIST = 0.0
      DO 20 I = 1,N
          DTLIST = DTLIST + STEP(I)
   20 CONTINUE
*
*       Set time-step list interval one-half the average step (N > 10).
      IF (N.GT.10) THEN
          DTLIST = 0.5*DTLIST/FLOAT(N)
      ELSE
          DTLIST = TCRIT
      END IF
      NNB = 1
   30 TLIST = TLIST + DTLIST
*
*       Make a list of particles due in the interval (0,TLIST).
      DO 40 J = 1,N
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              NNB = NNB + 1
              NLIST(NNB) = J
          END IF
   40 CONTINUE
*
*       Increase the interval for zero membership and save in NLIST(1).
      IF (NNB.EQ.1) GO TO 30
      NLIST(1) = NNB - 1
*
      RETURN
*
      END
