*
*             N B O D Y 1
*             ***********
*
*
*       Standard N-body code.
*       ---------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM NBODY1
*
      INCLUDE 'common1.h'
*
*
*       Initialize the timer.
      CALL CPUTIM(CPU0) 
*
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP
      IF (KSTART.EQ.1) GO TO 2
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
      CALL MYDUMP(0,1)
      IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine OUTPUT.
*
      CPU = TCOMP
      CPU0 = 0.0
      IF (KSTART.EQ.2) GO TO 4
*
*       Read modified parameters (KSTART = 3, 4 or 5).
      CALL MODIFY(KSTART)
      GO TO 4
*
    2 CPU = TCOMP
*
*       Read input parameters and perform initial setup.
      CALL START
*
*       Calculate total energy and produce output.
    3 CALL OUTPUT
*
*       Advance solutions until next output.
    4 CALL INTGRT
      GO TO 3
*
      END
