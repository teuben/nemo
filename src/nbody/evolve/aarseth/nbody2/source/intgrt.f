      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common2.h'
      COMMON/ISAVE/  LI0,LI,NS
      REAL*8  XI,YI,ZI
      REAL*4  XIDOT(3),TSLIST(NMAX+1)
*
*
*       Initialize first element of sorting array.
      TSLIST(1) = 0.0
*
*       Resume cycle after output (escape: IPHASE = -1; restart: LI0 = 0).
      IF (IPHASE.EQ.3.AND.LI0.GT.0) THEN
          IPHASE = 0
          NL = NLIST(1) + 1
          GO TO 1
      END IF
*
*       Form new time-step list after restart (#1 or 2), merge or escape. 
      IF (LI0.EQ.0.OR.IPHASE.EQ.-1) THEN
          IPHASE = 0
          TLIST = TIME
          GO TO 4
      END IF
*
*       Find next body to be advanced and set new time.
    1 LI = LI + 1
      IF (LI.GT.NL) GO TO 4
      I = NLIST(LI)  
      TIME = T0(I) + STEP(I)
*
*       Start integration cycle unless new time-step list is required.
      IF (TIME.LT.TLIST) GO TO 10
*
*       Form new time-step list and re-determine next body to be treated.
    4 NL = 1
      TLIST = TLIST + DTLIST
*
      DO 5 J = 1,N
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              NL = NL + 1
              TSLIST(NL) = T0(J) + STEP(J)
              NLIST(NL) = J 
          END IF
    5 CONTINUE
*
*       Increase interval and try again if no members found.
      IF (NL.EQ.1) THEN
          DTLIST = 2.0*DTLIST
          GO TO 4 
      END IF
      NLIST(1) = NL - 1 
*
*       Modify stabilizers by comparing insert membership with upper limit.
      IF (NS.LT.2) THEN
          STAB1 = MIN(1.1*STAB1,2.0)
          STAB2 = MIN(1.1*STAB2,4.0)
      ELSE
          STAB1 = MAX(0.9*STAB1,0.25)
          STAB2 = MAX(0.9*STAB2,0.5)
      END IF
*
*       Stabilize membership of NLIST in the range (STAB1,STAB2)*NNBMAX.
      IF (NL.GT.STAB2*NNBMAX) DTLIST = 0.75*DTLIST
      IF (NL.LT.STAB1*NNBMAX) DTLIST = 1.25*DTLIST
*
*       Sort new time-step list for sequential selection.
      CALL SORT2(NL,TSLIST,NLIST)
      LI = 1
      NS = 0
      GO TO 1
*
*       Advance the irregular step.
   10 CALL NBINT(I,IR,XI,YI,ZI,XIDOT)
*
      IF (IR.GT.0) THEN
*       Advance the regular step.
          CALL REGINT(I,XI,YI,ZI,XIDOT)
      END IF
*
*       Obtain new irregular step at end of integration cycle.
      CALL STEPI(I)
*
*       See whether current body is due before new NLIST loop.
      IF (TIME + STEP(I).LT.TLIST) THEN
*       Insert body #I in the correct sequential location.
          CALL INSERT(I,LI,NL,NS)
*       Ensure new NLIST if too many insertions.
          IF (NL.GT.NMAX - 10) THEN
              TLIST = TIME
          END IF
      END IF
*
*       See whether the merger indicator has been set.
      IF (IPHASE.GT.3) GO TO 100
*
*       Check next output time at the end of each integration cycle.
      IF (TIME.GT.TNEXT) THEN
          IPHASE = 3
*       Indicator for calling OUTPUT from MAIN.
          GO TO 100
      END IF
*
*       Advance counters and check timer & optional COMMON save.
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX) THEN
          NSTEPS = 0
          IF (KZ(1).GT.1) CALL MYDUMP(1,1)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          WRITE (6,30)
   30     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,40)  TIME, TCOMP, CPUTOT/60.0, ERRTOT
   40     FORMAT (//,9X,'COMMON SAVED AT TIME =',F8.2,3X,'TCOMP =',F8.2,
     &                            3X,'CPUTOT =',F7.1,3X,'ERRTOT =',F9.5)
      END IF
*
      STOP
*
  100 LI0 = LI
*
      RETURN
*
      END
