*
*             B R U T 4 -- GRAPE-free version
*             *******************************
*
*       Regularized N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 7.1.0 Oct 05).
*       -------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM NBODY4
*
      INCLUDE 'common4.h'
      EXTERNAL MERGE
*
*
*       Initialize the timer.
      CALL CPUTIM(CPU0)
*
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP, GPID
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
*         GPSTAT = 0
          CPU = TCOMP
          CALL START
          CALL ADJUST
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
*       Check restart counter (two extra restarts; set in routine CHECK).
          IF (NDUMP.GE.3) STOP
          CPU = TCOMP
          CPU0 = 0.0
*       Set IPHASE < -1 & ISEND < 0 for GRAPE initialization.
          IPHASE = -2
          ISEND = -1
*         GPSTAT = 0
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
          IF (KSTART.GT.2) THEN
              CALL MODIFY(KSTART)
          END IF
          CALL ZCNSTS(ZMET,ZPARS)
      END IF
*
*       Include optional movie initialization.
*     IF (KZ(40).GT.1) THEN
*         CALL MOVIE0
*     END IF
*
*       Advance solutions until next output or change of procedure.
    1 CALL INTGRT
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
          CALL KSREG
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
          CALL KSTERM
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          CALL ADJUST
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
          ISUB = 0
          CALL TRIPLE(ISUB)
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
          ISUB = 0
          CALL QUAD(ISUB)
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
          CALL MERGE
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
          CALL RESET
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
          ISUB = 0
          CALL CHAIN(ISUB)
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
