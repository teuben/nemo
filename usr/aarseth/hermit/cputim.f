      SUBROUTINE CPUTIM(TCOMP)
*
*
*       CPU time.
*       ---------
*
      REAL*8  TCOMP
      COMMON/ICPU0/ ICPU
      REAL*4  TARRAY(2)
*
*
*       Initialize timer (first call) or obtain elapsed time.
      IF (ICPU.EQ.0) THEN
*         CALL LIB$INIT_TIMER
          TCOMP = 0.0
*         TCOMP = ETIME(TARRAY)
*         TCOMP = MCLOCK()/6000.
          ICPU = 1
      ELSE
*         CALL LIB$STAT_TIMER(2,ITIME)
*         TCOMP = FLOAT(ITIME)/6000.0
          TCOMP = ETIME(TARRAY)/60.0
*         TCOMP = MCLOCK()/6000.
*       Elapsed CPU time in minutes on VAX, SUN or MIPS & IBM RS/6000.
      END IF
*
      RETURN
*
      END
