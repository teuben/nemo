*
*             H E R M I T 
*             ***********
*
*       N-body code with Hermite block-step integration.
*       ------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM HERMIT
*
      INCLUDE 'commonp.h'
*
*
*       Initialize the timer.
      CALL CPUTIM(CPU0)
*
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          CALL OUTPUT
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
          IF (KSTART.EQ.3) READ (5,*) J,K
          IF (J.GT.0) KZ(J) = K
      END IF
*
*       Advance solutions until next output or change of procedure.
    1 CALL INTGRT
*
*       Continue integration.
      GO TO 1
*
      END
