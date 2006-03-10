      SUBROUTINE SHRINK
*
*
*       Shrinking of large time-steps.
*       ------------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (DTM = 0.03125)
      REAL*8  RIJ(3),VIJ(3)
*
*
*       Determine impact to each high-velocity particle for large time-steps.
      DO 20 L = 1,NHI
          I = LISTV(L)
          DO 10 J = IFIRST,NTOT
*       Form relative quantities for large time-steps.
              IF (STEP(J).LT.DTM) GO TO 10
              RV = 0.0
              VV = 0.0
              DO 1 K = 1,3
                  RIJ(K) = X(K,I) - X(K,J)
                  VIJ(K) = XDOT(K,I) - XDOT(K,J)
                  RV = RV + RIJ(K)*VIJ(K)
                  VV = VV + VIJ(K)**2
    1         CONTINUE
*
*       Skip treatment for receding particles.
              IF (RV.GE.0.0D0) GO TO 10
*
*       Evaluate time of minimum approach and truncate to next time.
              DT = -RV/VV
              IT = 0
    2         DT = MIN(TNEXT(J) - TIME,DT)
*
*       Obtain minimum impact parameter for straight-line orbit.
              R2 = 0.0
              FJ2 = 0.0
              DO 5 K = 1,3
                  R2 = R2 + (RIJ(K) + VIJ(K)*DT)**2
                  FJ2 = FJ2 + F(K,J)**2
    5         CONTINUE
*
*       Compare force at minimum distance with total force on body #J (> 0).
              FI2 = (BODY(I)/R2)**2
              IF (FI2.LT.0.04*FJ2.OR.BODY(J).EQ.0.0D0) GO TO 10
*
*       See whether the time-step can be shortened.
              IF (T0(J) + 0.5*STEP(J).GT.TIME.AND.IT.LT.5) THEN
                  WRITE (29,8)  J, SQRT(R2), SQRT(FI2/FJ2), DT, STEP(J)
    8             FORMAT (' SHRINK    J R FI/FJ DT STEP ',
     &                                I6,F7.3,F6.2,1P,2E10.2)
                  CALL FLUSH(29)
                  STEP(J) = 0.5*STEP(J)
                  TNEXT(J) = T0(J) + STEP(J)
                  IF (IT.EQ.0) NSHORT = NSHORT + 1
                  IT = IT + 1
                  GO TO 2
              END IF
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
