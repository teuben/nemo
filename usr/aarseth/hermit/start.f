      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       ------------------------------------
*
      INCLUDE 'commonp.h'
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
*       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
*     CALL UNITS
*
*       Initialize sequential names & primary velocity.
      DO 20 I = 1,N
          NAME(I) = I
          DO 15 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Initialize fixed block steps (40 levels).
      CALL IBLOCK
*
*       Obtain force & first derivative for all bodies.
      CALL FPOLY1(1,N)
*
*       Set time-steps and initialize prediction variables.
      CALL STEPS(1,N)
*
      RETURN
*
      END
