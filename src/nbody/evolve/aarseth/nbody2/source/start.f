      SUBROUTINE START
*
*
*       Initialization procedures.
*       --------------------------
*
      INCLUDE 'common2.h'
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read and print input parameters.
      CALL INPUT
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
*
*       Check for optional creation of two separate subsystems.
      IF (KZ(17).GT.0) THEN
          CALL SUBSYS
      END IF
*
*       Set sequential name, maximum mass & total mass.
      BODY1 = 0.0
      ZMASS = 0.0
      DO 30 I = 1,N
          NAME(I) = I
          BODY1 = MAX(BODY1,BODY(I))
          ZMASS = ZMASS + BODY(I)
   30 CONTINUE
*
*       Define mean mass in scaled units.
      BODYM = ZMASS/FLOAT(N)
*
*       Initialize neighbour list & corresponding radius.
      RS0 = RC
      NNB = 0
      DO 40 I = 1,N
          CALL NBLIST(I,RS0)
          NNB = NNB + LIST(1,I)
   40 CONTINUE
*
*       Obtain force & first derivative for regular & irregular field.
      CALL FPOLY1(1,N)
*
*       Obtain second & third force derivatives and set time-steps.
      CALL FPOLY2(1,N)
*
*       Set new time-steps and initialize divided differences.
      CALL STEPS(1,N)
*
*       Initialize the time-step list used to find next body.
      DTLIST = 100.0
      DO 50 I = 1,N
          DTLIST = MIN(DTLIST,STEP(I))
   50 CONTINUE
*
*       Set initial time-step list interval twice the smallest step.
      DTLIST = 2.0*DTLIST
*
*       Check the average neighbour number.
      ZNB = FLOAT(NNB)/FLOAT(N)
      IF (ZNB.LT.0.25*ZNBMAX) THEN
          WRITE (6,60)  ZNB
   60     FORMAT (/,12X,'WARNING!   SMALL NEIGHBOUR NUMBERS   <NNB> =',
     &                                                             F5.1)
      END IF
*
      RETURN
*
      END
