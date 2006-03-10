      SUBROUTINE KEPLER(I,DTI)
*
*
*       Step reduction of hierarchy.
*       ----------------------------
*
*
      INCLUDE 'common4.h'
      COMMON/EXTRA/  RPERT2,I1,I2
*
*
*       Only examine hierarchical configurations (NP < 5 OR STEP < DTMIN).
      IF (LIST(1,I1).GT.4.AND.DTI.GT.DTMIN) GO TO 30
*
*       Set perturber membership & Kepler period.
      NP1 = LIST(1,I1) + 1
      SEMI = -0.5D0*BODY(I)/H(I-N)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Consider the most dominant perturbers having comparable step to c.m.
      DO 20 L = 2,NP1
          J = LIST(L,I1)
          IF (STEP(J).GT.3.0*STEP(I)) GO TO 20 
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          RIJ = SQRT(RIJ2)
          DT2 = 0.14*RIJ*SQRT(ETA*RIJ/(BODY(I) + BODY(J)))
          DT = 0.25D0*SQRT(ETA)*RIJ*TK/SEMI
*       Compare predicted c.m. step with conservative Kepler expressions.
          DT = MIN(DT2,DT)
          IF (DTI.LT.DT) GO TO 20 
          DTI = DT
*
*       Check whether to reduce step of dominant perturber.
          IF (STEP(J).LT.DT) GO TO 20
          IF (KZ(36).LE.0) THEN
              A2 = T0(J) + STEP(J)
*       Adopt maximum reduction factor 0.5 ensuring T0 + STEP > TIME.
              STEP(J) = STEP(J) - 0.5*(A2 - TIME)
          ELSE IF (T0(J).EQ.TIME) THEN
              STEP(J) = 0.5D0*STEP(J)
              TNEXT(J) = T0(J) + STEP(J)
              GO TO 20
          END IF
*
*       See whether body #J should be added to NLIST.
          IF (T0(J) + STEP(J).LT.TLIST) THEN
              CALL NLMOD(J,1)
          END IF
   20 CONTINUE
*
   30 RETURN
*
      END
