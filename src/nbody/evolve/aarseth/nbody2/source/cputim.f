      SUBROUTINE CPUTIM(TCOMP)
*
*
*       CPU time.
*       ---------
*
      COMMON/ICPU0/ ICPU
      REAL*4  TARRAY(2)
*
*
*       Initialize timer (first call) or obtain elapsed time.
      IF (ICPU.EQ.0) THEN
*         CALL LIB$INIT_TIMER
*         TCOMP = 0.0
          TCOMP = ETIME(TARRAY)
*         TCOMP = MCLOCK()/6000.0
          ICPU = 1
      ELSE
*         CALL LIB$STAT_TIMER(2,ITIME)
*         TCOMP = FLOAT(ITIME)/6000.0
*       Elapsed CPU time in minutes on VAX.
          TCOMP = ETIME(TARRAY)/60.0
*       Elapsed CPU time in minutes on SUN or MIPS.
*         TCOMP = MCLOCK()/6000.0
*       Elapsed CPU time in minutes on IBM workstation.
      END IF
*
      RETURN
*
      END
**
** uncomment this if you need a dummy ETIME; not all compilers have this
**
      REAL FUNCTION ETIME(TARRAY)
      REAL TARRAY(2)
      TARRAY(1) = 1.0
      TARRAY(2) = 2.0
      ETIME = 1.0
      RETURN
      END

