      SUBROUTINE LEVELS
*
*
*       Diagnostics for block time-step levels.
*       ---------------------------------------
*
      INCLUDE 'common4.h'
      INTEGER  IHIST(40)
*
*
*       Initialize histogram counter.
      DO 10 J = 1,40
          IHIST(J) = 0
   10 CONTINUE
*
*       Loop over all single particles & c.m.
      JMAX = 0
      FAC = 1.0/LOG(1.9999999)
      DO 20 I = IFIRST,NTOT
          IF (BODY(I).EQ.0.0D0) GO TO 20
          DT = STEP(I)
*       Adopt fast procedure for first four levels.
          IF (DT.GE.0.125) THEN
              IF (DT.GE.0.5) THEN
                  IF (DT.EQ.1.0) THEN
                      J = 1
                   ELSE
                       J = 2
                   END IF
              ELSE IF (DT.EQ.0.25) THEN
                  J = 3
              ELSE
                  J = 4
              END IF
*       Consider the next four levels in a similar way.
          ELSE IF (DT.GE.0.0078125) THEN
              IF (DT.GE.0.03125) THEN
                  IF (DT.EQ.0.0625) THEN
                      J = 5
                  ELSE
                      J = 6
                  END IF
              ELSE IF (DT.EQ.0.015625) THEN
                  J = 7
              ELSE
                  J = 8
              END IF
*      Include level 9 & 10 explicitly.
          ELSE IF (DT.GE.0.001953125) THEN
              IF (DT.EQ.0.00390625) THEN
                  J = 9
              ELSE
                  J = 10
              END IF
          ELSE
*      Treat the tail of step distribution by actual LOG.
          J = 1 - LOG(DT)*FAC
          END IF
          IHIST(J) = IHIST(J) + 1
          JMAX = MAX(J,JMAX)
   20 CONTINUE
*
*       Print histogram of block-steps.
      WRITE (6,30)  (IHIST(J),J=1,JMAX)
   30 FORMAT (' #6  STEP   ',10I6,/,12X,15I6)
*
*       Check optional diagnostics of active pipes.
      IF (KZ(33).GT.1) THEN
          WRITE (6,40)  (IPIPE(J),J=1,13)
   40     FORMAT (' #7  IPIPE  ',2I11,2I10,2I9,7I8)
*       Rectify any counters > 2^{32}.
          DO 50 J = 1,20
              IF (IPIPE(J).GT.2000000000.OR.IPIPE(J).LT.0) THEN
                  IPIPE(J) = 0
              END IF
   50     CONTINUE
      END IF
*
      RETURN
*
      END
