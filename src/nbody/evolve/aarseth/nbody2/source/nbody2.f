*
*             N B O D Y 2
*             ***********
*
*       Ahmad-Cohen N-body code.
*       ------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM NBODY2
*
      INCLUDE 'common2.h'
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
*       Safety indicator preventing repeated restarts set in routine OUTPUT.
*
          CPU = TCOMP
          CPU0 = 0.0
*       Set IPHASE = -1 to ensure new NLIST in routine INTGRT.
          IPHASE = -1
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
      END IF
*
*       Advance solutions until next output or change of procedure.
    1 CALL INTGRT
*
*       Check merger option & indicator.
      IF (KZ(12).GT.0.AND.IPHASE.EQ.4) THEN
          CALL MERGE
          GO TO 1
      END IF
*
*       Calculate total energy and produce main output.
      CALL OUTPUT
*
      GO TO 1 
*
      END
