      SUBROUTINE INSERT(I,LI,NL,NS)
*
*
*       Insert particle index in time-step list.
*       ----------------------------------------
*
      INCLUDE 'common2.h'
      REAL*8  TI,TI1,TI2,TJ
*
*
*       Avoid increasing NLIST if body #I can be exchanged with next member.
      LI2 = MIN(LI + 2,NL)
      I2 = NLIST(LI2)
      TI = TIME + STEP(I)
      TI2 = T0(I2) + STEP(I2)
*
*       Check swapping condition (#I due before #I2 but after #I1 at LI + 1).
      IF (TI.LE.TI2) THEN
          LI1 = MIN(LI + 1,NL)
          I1 = NLIST(LI1)
          TI1 = T0(I1) + STEP(I1)
          IF (TI.GT.TI1) THEN
              NLIST(LI) = NLIST(LI1)
              NLIST(LI1) = I
          END IF
*       Reduce pointer so that current location will be selected again.
          LI = LI - 1
          GO TO 50
*       Also swap if body #I is among the two last members and TI > TI2.
      ELSE IF (LI.GT.NL - 2) THEN
          NLIST(LI) = NLIST(LI2)
          NLIST(LI2) = I
          LI = LI - 1
          GO TO 50
      END IF
*
*       Estimate the insert index from the remaining time interval.
      FAC = STEP(I)/(TLIST - TIME)
      LSTAR = LI + FLOAT(NL - LI)*FAC
*
*       Improve insert index by another iteration (check LI < LSTAR <= NL).
      J = NLIST(LSTAR)
      TJ = T0(J) + STEP(J)
      FAC = (TI - TJ)/(TLIST - TJ)
      LSTAR = LSTAR + FLOAT(NL - LSTAR)*FAC
      LSTAR = MAX(LSTAR,LI + 1)
      LSTAR = MIN(LSTAR,NL)
      J = NLIST(LSTAR)
*
*       Determine correct index by comparing neighbouring NLIST members.
      IF (TI.GE.T0(J) + STEP(J)) THEN
          L1 = LSTAR + 1
          LSTAR = L1
          DO 10 L = L1,NL
              J = NLIST(L)
              TJ = T0(J) + STEP(J)
*       Advance index until TI < TJ or last member.
              IF (TI.LT.TJ) GO TO 30 
              LSTAR = L + 1
   10     CONTINUE
      ELSE
   20     LSTAR = LSTAR - 1
          J = NLIST(LSTAR)
          IF (TI.LT.T0(J) + STEP(J)) GO TO 20
          LSTAR = LSTAR + 1
*       Ensure that current index LI is not chosen due to round-off.
          LSTAR = MAX(LSTAR,LI + 1)
      END IF 
*
*       Create free location at LSTAR by moving all subsequent members down.
   30 DO 40 L = NL,LSTAR,-1
          NLIST(L+1) = NLIST(L)
   40 CONTINUE 
*
*       Insert body #I and update memberships.
      NLIST(LSTAR) = I
      NLIST(1) = NL
      NL = NL + 1
      NS = NS + 1
*
   50 RETURN
*
      END
